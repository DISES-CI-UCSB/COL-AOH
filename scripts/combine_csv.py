import os
import pandas as pd
import glob

def combine_csv_files(root_folder):
    # Use glob to get all CSV files in root_folder and its subfolders
    all_files = glob.glob(os.path.join(root_folder, "**/*.csv"), recursive=True)
    
    # Create an empty list to store dataframes
    df_list = []
    
    # Read each CSV file and append to the list
    for file in all_files:
        try:
            df = pd.read_csv(file)
            # Optionally add a column to identify the source file
            df['source_file'] = os.path.basename(file)
            df_list.append(df)
            print(f"Successfully read: {file}")
        except Exception as e:
            print(f"Error reading {file}: {str(e)}")
    
    if not df_list:
        print("No CSV files found in the specified directory and its subfolders.")
        return
    
    # Combine all dataframes
    combined_df = pd.concat(df_list, ignore_index=True)
    
    # Create output filename
    output_file = os.path.join(root_folder, "combined_output.csv")
    
    # Save the combined dataframe
    combined_df.to_csv(output_file, index=False)
    print(f"\nCombined CSV has been saved to: {output_file}")
    print(f"Total number of files combined: {len(df_list)}")

if __name__ == "__main__":
    # Get the folder path from user input
    folder_path = input("Enter the folder path containing CSV files: ")
    
    # Check if the folder exists
    if os.path.exists(folder_path):
        combine_csv_files(folder_path)
    else:
        print("The specified folder does not exist. Please check the path and try again.") 