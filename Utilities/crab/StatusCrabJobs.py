import os
import subprocess
import argparse

def run_crab_status(base_dir):
    # Convert base directory to absolute path
    base_dir = os.path.abspath(base_dir)
    print(f"Base directory: {base_dir}")  # Debugging print
    
    # Check if the base directory exists
    if not os.path.isdir(base_dir):
        print(f"Error: Directory '{base_dir}' does not exist.")
        return
    
    # Iterate over all items in the base directory
    found_dirs = False
    for item in os.listdir(base_dir):
        sub_dir = os.path.join(base_dir, item)
        
        # Check if the item is a directory and starts with "Crab"
        if os.path.isdir(sub_dir) and item.startswith("crab"):
            found_dirs = True
            relative_path = os.path.relpath(sub_dir, os.getcwd())
            print(f"Found matching directory: {relative_path}")  # Debugging print
            print(f"Running 'crab status --verboseErrors -d {relative_path}'")
            try:
                # Run the command and capture the output
                result = subprocess.run(["crab", "status","--verboseErrors", "-d", relative_path], capture_output=True, text=True)
                print(result.stdout)
                print(result.stderr)
            except Exception as e:
                print(f"Failed to run command for {relative_path}: {e}")

    if not found_dirs:
        print("No subdirectories starting with 'Crab' were found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run 'crab status' for all subdirectories in the specified folder.")
    parser.add_argument("base_folder", help="Path to the base folder containing subdirectories.")
    args = parser.parse_args()
    
    run_crab_status(args.base_folder)
