# Load necessary libraries
library(sf)
library(terra)
library(tools) # For file_path_sans_ext
library(dplyr) # For data manipulation, especially mutate and select

# --- Main Function to Perform All Spatial Appends ---
# attributes_to_transfer is now a character vector
perform_all_spatial_appends <- function(point_layer_path, raster_folders, polygon_folders, attributes_to_transfer) {
  
  # --- 1. Setup Paths and Load Initial Point Layer ---
  message(paste0("Starting processing for: ", point_layer_path))
  
  # Dynamically create output paths
  input_path_dir <- dirname(point_layer_path)
  input_file_name <- file_path_sans_ext(basename(point_layer_path))
  
  # Output path for CSV
  output_csv_path <- file.path(input_path_dir, paste0(input_file_name, "_output.csv"))
  # Output path for GeoPackage
  output_gpkg_path <- file.path(input_path_dir, paste0(input_file_name, "_output.gpkg"))
  
  # Load point layer
  point_layer <- st_read(point_layer_path, quiet = TRUE)
  if (inherits(point_layer, "try-error") || nrow(point_layer) == 0) {
    stop("Error: Point layer is invalid or empty.")
  }
  
  # --- IMPORTANT CHANGE: Drop M (and Z if desired) dimensions from point layer ---
  # st_zm(drop = TRUE) will drop the M dimension. If it's XYZM, it becomes XYZ.
  # If you want to force 2D (XY) even if Z is present, use st_zm(drop = TRUE, what = "ZM").
  # For point features and spatial joins, XY is usually sufficient.
  point_layer <- st_zm(point_layer, drop = TRUE) 
  
  point_crs <- st_crs(point_layer)
  message(paste0("Point layer loaded. CRS: ", point_crs$input))
  
  # Initialize the result layer with the original point layer
  result_layer <- point_layer
  
  # --- 2. Process Raster Data ---
  message("\n--- Starting Raster Data Processing ---")

  if (length(raster_folders) == 0 || all(sapply(raster_folders, function(f) !dir.exists(f)))) {
    message("No valid raster folders provided. Skipping raster processing.")
  } else {
    for (current_raster_folder in raster_folders) {
      if (!dir.exists(current_raster_folder)) {
        warning(paste0("Raster folder does not exist: ", current_raster_folder, ". Skipping."))
        next
      }
      folder_base_name <- basename(current_raster_folder) # Get the folder name

      raster_files <- list.files(current_raster_folder, pattern = "\\.tif$|\\.asc$", full.names = TRUE)

      if (length(raster_files) == 0) {
        message(paste0("No raster files found in folder: ", current_raster_folder, ". Skipping."))
        next
      }

      for (raster_path in raster_files) {
        raster_file_name <- file_path_sans_ext(basename(raster_path))
        # New prefix format: folder_name_file_name_
        column_prefix <- paste0(folder_base_name, "_", raster_file_name, "_")
        message(paste0("Processing raster: ", raster_file_name, " from folder: ", folder_base_name))

        raster_to_sample <- rast(raster_path)
        if (is.null(raster_to_sample)) {
          warning(paste0("Failed to load raster: ", raster_path, ". Skipping."))
          next
        }

        # --- IMPORTANT CHANGE FOR MISSING RASTER CRS ---
        raster_crs_info <- crs(raster_to_sample, describe = TRUE) # Get detailed CRS info

        if (is.na(raster_crs_info$code) || raster_crs_info$code == "" || raster_crs_info$code == "undef") { # Check for 'undef' as well
          message(paste0("  Raster '", raster_file_name, "' has no CRS defined. Assuming it's in the same CRS as the point layer (", point_crs$input, ")."))
          # Assign the point layer's CRS to the raster
          crs(raster_to_sample) <- point_crs$wkt
          # Now, the raster has a CRS assigned, so the next comparison/reprojection will work.
        }

        # Reproject raster if CRS differs from point layer (now including cases where it was just assigned)
        if (crs(raster_to_sample) != point_crs$wkt) {
          message(paste0("  Reprojecting raster '", raster_file_name, "' to match point layer CRS."))
          raster_to_sample <- project(raster_to_sample, point_crs$wkt, method = "near")
        }

        # Extract values
        extracted_values <- terra::extract(raster_to_sample, result_layer)

        if (ncol(extracted_values) > 1) {
          # Make sure 'ID' column is for merging, not to be prefixed
          names_to_prefix <- setdiff(names(extracted_values), "ID")

          # Handle cases where terra::extract might return band names like "X1", "X2", etc.
          # We want to use the actual band names if available, or just append the prefix if they are generic.
          if(length(names(raster_to_sample)) == length(names_to_prefix)) { # Check if extracted names match original band names
            original_band_names <- names(raster_to_sample)
          } else { # Fallback if names from extract are not easily mapped to original band names
            original_band_names <- names_to_prefix # Use the extracted names directly
          }

          # Rename columns with the new prefix
          renamed_extracted_values <- extracted_values %>%
            select(ID, all_of(names_to_prefix)) %>% # Ensure 'ID' is kept
            rename_with(~ paste0(column_prefix, original_band_names[match(.x, names_to_prefix)]), .cols = all_of(names_to_prefix))

          # Ensure result_layer has a temporary ID for merging
          if (!"temp_id_for_merge" %in% names(result_layer)) {
            result_layer$temp_id_for_merge <- 1:nrow(result_layer)
          }
          renamed_extracted_values$temp_id_for_merge <- renamed_extracted_values$ID

          # Perform a non-spatial join to add the new raster attributes
          result_layer <- result_layer %>%
            left_join(st_drop_geometry(renamed_extracted_values %>% select(-ID)), # Exclude terra's ID
                      by = "temp_id_for_merge")

        } else {
          # If no values were extracted (e.g., points outside raster extent)
          message(paste0("No raster values extracted for ", raster_file_name, ". Adding NA column(s)."))
          # Add NA column(s) based on expected raster layers, using the new prefix
          expected_cols <- paste0(column_prefix, names(raster_to_sample))
          for (col_name in expected_cols) {
            if (!col_name %in% names(result_layer)) {
              result_layer[[col_name]] <- NA
            }
          }
        }
        rm(raster_to_sample, extracted_values, renamed_extracted_values) # Clean up to free memory
      }
    }
    # Remove the temporary merge ID after all raster joins are done
    if ("temp_id_for_merge" %in% names(result_layer)) {
      result_layer <- result_layer %>% select(-temp_id_for_merge)
    }
  }
  
  # --- 3. Process Vector Data (Spatial Join) ---
  message("\n--- Starting Vector Data Spatial Join ---")
  
  if (length(polygon_folders) == 0 || all(sapply(polygon_folders, function(f) !dir.exists(f)))) {
    message("No valid polygon folders provided. Skipping vector processing.")
  } else {
    # Ensure temp_id_for_merge_v is present ONCE for the whole vector processing
    # This ID will be used to correctly join back attributes after each spatial join.
    if (!"temp_id_for_merge_v" %in% names(result_layer)) {
      result_layer$temp_id_for_merge_v <- 1:nrow(result_layer)
    }
    
    for (current_polygon_folder in polygon_folders) {
      if (!dir.exists(current_polygon_folder)) {
        warning(paste0("Polygon folder does not exist: ", current_polygon_folder, ". Skipping."))
        next
      }
      folder_base_name <- basename(current_polygon_folder) # Get the folder name
      
      # Adjusted to also look for .gpkg files for polygons
      polygon_files <- list.files(current_polygon_folder, pattern = "\\.shp$|\\.gpkg$", full.names = TRUE)
      
      if (length(polygon_files) == 0) {
        message(paste0("No files found in folder: ", current_polygon_folder, ". Skipping."))
        next
      }
      
      for (polygon_path in polygon_files) {
        polygon_file_name <- basename(polygon_path)
        polygon_base_name <- file_path_sans_ext(polygon_file_name)
        message(paste0("Processing polygon layer: ", polygon_file_name, " from folder: ", folder_base_name))
        
        polygon_layer <- st_read(polygon_path, quiet = TRUE)
        if (inherits(polygon_layer, "try-error") || nrow(polygon_layer) == 0) {
          warning(paste0("Error: Could not load or empty polygon layer: ", polygon_file_name, ". Skipping."))
          rm(polygon_layer) # Ensure cleanup
          next
        }
        
        # --- IMPORTANT CHANGE: Drop M (and Z if desired) dimensions from polygon layer ---
        polygon_layer <- st_zm(polygon_layer, drop = TRUE) 
        
        # Ensure CRS matches
        if (st_crs(polygon_layer) != point_crs) {
          message(paste0("Reprojecting polygon layer '", polygon_file_name, "' to match point layer CRS."))
          polygon_layer <- st_transform(polygon_layer, point_crs)
        }
        
        # Select only the attributes to transfer and the geometry from the polygon layer
        # Using `any_of()` here is safer if some attributes might be missing in some shapefiles
        polygon_for_join <- polygon_layer %>% select(any_of(attributes_to_transfer))
        
        # Check if any of the requested attributes were actually found in this polygon layer
        found_attrs_in_polygon <- intersect(attributes_to_transfer, names(st_drop_geometry(polygon_for_join)))
        if (length(found_attrs_in_polygon) == 0) {
          warning(paste0("None of the specified attributes (", paste(attributes_to_transfer, collapse = ", "), ") found in polygon layer: ", polygon_file_name, ". Skipping join for this layer."))
          rm(polygon_layer, polygon_for_join)
          next # Skip to next polygon file
        }
        
        # Perform the spatial join. `st_join` directly on `result_layer` will keep its geometry and attributes.
        # It handles column name conflicts by appending .y
        temp_joined_layer <- st_join(result_layer, polygon_for_join, left = TRUE)
        
        # Identify the columns that were newly added by this join from polygon_for_join
        # These will be the names from found_attrs_in_polygon, potentially with '.y' suffix
        newly_added_col_names <- setdiff(names(temp_joined_layer), names(result_layer)) # Find truly new columns
        
        # Filter to only the ones that match our expected transferred attributes (considering .y suffix)
        cols_to_rename_actual <- c()
        for(attr_name in found_attrs_in_polygon) {
          current_name_in_joined <- NA
          if (attr_name %in% newly_added_col_names) { # Case: name didn't clash and is new
            current_name_in_joined <- attr_name
          } else if (paste0(attr_name, ".y") %in% newly_added_col_names) { # Case: name clashed
            current_name_in_joined <- paste0(attr_name, ".y")
          }
          if (!is.na(current_name_in_joined)) {
            cols_to_rename_actual <- c(cols_to_rename_actual, current_name_in_joined)
          }
        }
        
        # Rename these columns with the desired prefix
        if (length(cols_to_rename_actual) > 0) {
          # Create a dataframe with just the temp ID and the columns to rename
          # We use st_drop_geometry to work with data frame for joining
          attrs_to_update <- st_drop_geometry(temp_joined_layer) %>%
            select(temp_id_for_merge_v, all_of(cols_to_rename_actual)) %>%
            rename_with(~ paste0(folder_base_name, "_", polygon_base_name, "_", 
                                 sub("\\.y$", "", .x)), # Remove '.y' if present before adding prefix
                        .cols = all_of(cols_to_rename_actual))
          
          # Merge these renamed attributes back to the result_layer based on temp_id_for_merge_v
          # This ensures only the *newly prefixed* columns are added, without duplicating original ones.
          result_layer <- result_layer %>%
            left_join(attrs_to_update, by = "temp_id_for_merge_v")
          
        } else {
          # If no columns were added or matched (e.g., all points outside polygons),
          # explicitly add NA columns for requested attributes to ensure consistency.
          message(paste0("  No values joined for ", polygon_file_name, ". Adding NA columns for requested attributes."))
          for (attr in attributes_to_transfer) { # Add for all requested attributes
            new_col_name <- paste0(folder_base_name, "_", polygon_base_name, "_", attr)
            if (!new_col_name %in% names(result_layer)) { # Only add if not already present from previous join
              result_layer[[new_col_name]] <- NA
            }
          }
        }
        
        # Clean up
        rm(polygon_layer, polygon_for_join, temp_joined_layer, attrs_to_update)
      }
    }
    # Remove the temporary merge ID after all vector joins are done
    if ("temp_id_for_merge_v" %in% names(result_layer)) {
      result_layer <- result_layer %>% select(-temp_id_for_merge_v)
    }
  }
  
  # --- 4. Save Final Results as CSV and GeoPackage ---
  message("\n--- Saving Final Outputs ---")
  if (!is.null(result_layer) && nrow(result_layer) > 0) {
    # Save as CSV
    coords <- as.data.frame(st_coordinates(result_layer))
    colnames(coords) <- c("X_coordinate", "Y_coordinate")
    result_df <- st_drop_geometry(result_layer)
    final_output_df <- cbind(coords, result_df)
    write.csv(final_output_df, output_csv_path, row.names = FALSE)
    message(paste0("  Final CSV saved to: ", output_csv_path))
    
    # Save as GeoPackage
    st_write(result_layer, output_gpkg_path, driver = "GPKG", delete_layer = TRUE, append = FALSE, quiet = TRUE)
    message(paste0("  Final GeoPackage saved to: ", output_gpkg_path))
    
  } else {
    stop("Error: Final result layer is empty or invalid. No output saved.")
  }
  
  return(result_layer) # Return the sf object for potential further R use
}

# --- Define Input Parameters for Both Cases ---

# Case 1: Fassadenpunkte
point_layer_path_1 <- 'D:/Gesundheitsdaten_3_R/Daten/Mainz-Bingen/Fassadenpunkte_corrected.gpkg'
raster_folders_1 <- c('D:/Gesundheitsdaten_3_R/Daten/Fluglärm/', 'D:/Gesundheitsdaten_3_R/Daten/OtherRasterData/') # Can be one or multiple folders
polygon_folders_1 <- c('D:/Gesundheitsdaten_3_R/Daten/Klimprax/', 'D:/Gesundheitsdaten_3_R/Daten/OtherVectorData/') # Can be one or multiple folders
# Now attributes_to_transfer is a VECTOR of column names
attributes_to_transfer_1 <- c("GRID_CODE", "GEN", "plz_code")  

# Case 2: Strasse
point_layer_path_2 <- 'D:/Gesundheitsdaten_3_R/Daten/Mainz/strasse_2022_corrected.gpkg'
raster_folders_2 <- c('D:/Gesundheitsdaten_3_R/Daten/Fluglärm/', 'D:/Gesundheitsdaten_3_R/Daten/OtherRasterData/') # Can be one or multiple folders
polygon_folders_2 <- c('D:/Gesundheitsdaten_3_R/Daten/Klimprax/', 'D:/Gesundheitsdaten_3_R/Daten/OtherVectorData/', 'D:/Gesundheitsdaten_3_R/Daten/Mainz_VectorData/') # Can be one or multiple folders
# Now attributes_to_transfer is a VECTOR of column names
attributes_to_transfer_2 <- c("GRID_CODE", "GEN", "plz_code", "LP1", "LP2", "LP3", "LP4")  


# --- Execute the Combined Function ---
message("\n### Processing Fassadenpunkte ###")
final_fassadenpunkte_layer <- perform_all_spatial_appends(
  point_layer_path = point_layer_path_1,
  raster_folders = raster_folders_1,
  polygon_folders = polygon_folders_1,
  attributes_to_transfer = attributes_to_transfer_1 # Pass the vector of attributes
)

message("\n### Processing Strasse ###")
final_strasse_layer <- perform_all_spatial_appends(
  point_layer_path = point_layer_path_2,
  raster_folders = raster_folders_2,
  polygon_folders = polygon_folders_2,
  attributes_to_transfer = attributes_to_transfer_2 # Pass the vector of attributes
)

# You can still use 'final_fassadenpunkte_layer' and 'final_strasse_layer'
# as sf objects in R for further spatial analysis if needed,
# even though the primary output is now CSV.