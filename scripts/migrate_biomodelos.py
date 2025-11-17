import os
import glob
import shutil


def copy_maxent_files(source_dir, destination_dir):
    # Create destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Counter for copied files
    copied_count = 0

    subfolders = ['peces']

    for folder in subfolders:
        dest_folder = os.path.join(destination_dir, folder)
        if not os.path.exists(dest_folder):
            os.makedirs(dest_folder)

        folder = os.path.join(source_dir, folder)
        subdirs = [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]
        for thisdir in subdirs:
            path_dir = os.path.join(folder, thisdir)
            species_name = thisdir.replace(".", "_")
            maxent_file_path = os.path.join(path_dir, 'ensembles/current/MAXENT', species_name + "_MAXENT.tif")

            if os.path.exists(maxent_file_path):
                dest_path = os.path.join(dest_folder, species_name + ".tif")

                # Copy the file
                try:
                    shutil.copy2(maxent_file_path, dest_path)
                    print(f"Successfully copied to: {dest_path}")
                    copied_count += 1
                except Exception as e:
                    print(f"Error copying {file}: {str(e)}")

    print(f"\nTotal files copied: {copied_count}")

def copy_occ_files(source_dir, destination_dir):
    # Create destination directory if it doesn't exist
    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    # Counter for copied files
    copied_count = 0

    subfolders = ['anfibios', 'aves', 'mamiferos', 'squamata']

    for folder in subfolders:
        dest_folder = os.path.join(destination_dir, folder)
        if not os.path.exists(dest_folder):
            os.makedirs(dest_folder)

        folder = os.path.join(source_dir, folder)
        subdirs = [d for d in os.listdir(folder) if os.path.isdir(os.path.join(folder, d))]
        for thisdir in subdirs:
            path_dir = os.path.join(folder, thisdir)
            species_name = thisdir.replace(".", "_")
            occ_file_path = os.path.join(path_dir, 'occurrences/formated_occ.csv')

            if os.path.exists(occ_file_path):
                dest_path = os.path.join(dest_folder, species_name + ".csv")

                # Copy the file
                try:
                    shutil.copy2(occ_file_path, dest_path)
                    print(f"Successfully copied to: {dest_path}")
                    copied_count += 1
                except Exception as e:
                    print(f"Error copying {file}: {str(e)}")

    print(f"\nTotal files copied: {copied_count}")


# Example usage
source_directory = "/home/scale-lab/dises/biodiversity/biomodelos/full_models"
destination_directory = "/home/wenxinyang/Desktop/wenxinyang/colander/biomodelos/occ_pts"

if __name__ == "__main__":
    # copy_maxent_files(source_directory, destination_directory)
    copy_occ_files(source_directory, destination_directory)