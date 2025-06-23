import matplotlib.pyplot as plt

def parse_bvh(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()

    joint_channels = []
    motion_start = None
    channel_names = []
    joint_name_stack = []

    # === Step 1: Parse HIERARCHY ===
    for i, line in enumerate(lines):
        stripped = line.strip()

        if stripped.startswith("MOTION"):
            motion_start = i
            break

        if stripped.startswith("ROOT") or stripped.startswith("JOINT"):
            tokens = stripped.split()
            joint_name_stack.append(tokens[1])
        elif stripped.startswith("End Site"):
            joint_name_stack.append("EndSite")
        elif stripped.startswith("CHANNELS"):
            tokens = stripped.split()
            num_channels = int(tokens[1])
            joint_name = joint_name_stack[-1]
            channels = tokens[2:2+num_channels]
            joint_channels.append((joint_name, channels))
            channel_names += [(joint_name, ch) for ch in channels]
        elif stripped.startswith("}"):
            if joint_name_stack:
                joint_name_stack.pop()

    # === Step 2: Parse MOTION Header ===
    num_frames = int(lines[motion_start+1].strip().split()[-1])
    frame_time = float(lines[motion_start+2].strip().split()[-1])
    data_lines = lines[motion_start+3:]

    assert len(data_lines) == num_frames

    motion_data = []
    for line in data_lines:
        motion_data.append([float(x) for x in line.strip().split()])

    return channel_names, motion_data, frame_time

# === Step 3: 提取某个关节的角度 ===
def extract_joint_rotation(channel_names, motion_data, joint_name, rotation_type="Xrotation"):
    for i, (jname, cname) in enumerate(channel_names):
        if jname == joint_name and cname == rotation_type:
            index = i
            break
    else:
        raise ValueError(f"Cannot find {rotation_type} for joint {joint_name}")

    return [frame[index] for frame in motion_data]

# === Step 4: 使用示例 ===
bvh_file_path = "/home/joe/Desktop/MASS (copy)/data/motion/walk.bvh"
channel_names, motion_data, frame_time = parse_bvh(bvh_file_path)

# ✅ 同时提取左右膝盖的 Xrotation
right_knee = extract_joint_rotation(channel_names, motion_data, "Character1_RightLeg", "Xrotation")
left_knee = extract_joint_rotation(channel_names, motion_data, "Character1_LeftLeg", "Xrotation")

# === Step 5: 绘图 ===
plt.figure(figsize=(10, 4))
plt.plot(right_knee, label="Right Knee Xrotation", color='blue')
plt.plot(left_knee, label="Left Knee Xrotation", color='red')
plt.xlabel("Frame")
plt.ylabel("Angle (deg)")
plt.title("Left & Right Knee Angle Trajectory from BVH")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()
