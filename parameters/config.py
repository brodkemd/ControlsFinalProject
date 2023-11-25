import os, json

root = os.path.dirname(os.path.dirname(__file__))
with open(os.path.join(root, "data", "user_config.json")) as f:
    data = json.loads(f.read())

matlab = data["matlab"]

if not os.path.exists(matlab):
    print("Matlab path in config does not exist:", matlab)
    exit()

del data
print(root, matlab)