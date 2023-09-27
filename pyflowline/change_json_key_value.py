import json

def change_json_key_value(sFilename_json_in, sKey, new_value):
    # Read the JSON file
    with open(sFilename_json_in, 'r') as file:
        data = json.load(file)

    # Update the value associated with the specified sKey
    data[sKey] = new_value

    # Write the updated data back to the JSON file
    with open(sFilename_json_in, 'w') as file:
        json.dump(data, file, indent=4)