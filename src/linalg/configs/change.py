import os

def update_json_files():
    current_directory = os.getcwd()
    json_files = [f for f in os.listdir(current_directory) if f.endswith('.json')]

    for json_file in json_files:
        file_path = os.path.join(current_directory, json_file)
        update_json_file(file_path)
        print(f"Updated values in {json_file}")

def update_json_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    for i in range(len(lines)):
        if "print" in lines[i]:
            # Assuming the line is long enough to have a third-to-last character
            if len(lines[i]) > 4:
                lines[i] = lines[i][:-4] + '0' + lines[i][-3:]

    with open(file_path, 'w') as file:
        file.writelines(lines)

if __name__ == "__main__":
    update_json_files()
