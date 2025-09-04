import os

objects = [
    "sertraline_exh32_r2",
    "indatraline_exh32_r1",
    "leelamine_exh32_r2",
]

for obj_name in objects:
    base_path = rf"C:\Users\lam16\Documents\Cours\2025\projet-long\martini_to_allatom\docking\vina\results\{obj_name}"
    base_filename = f"{obj_name}_fr"

    num_states = cmd.count_states(obj_name)

    for state in range(1, num_states + 1):
        filename = os.path.join(base_path, f"{base_filename}{state}.pdb")
        cmd.save(filename, obj_name, state=state)
        print(f"Saved frame {state} as {filename}")
