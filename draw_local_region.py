import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def read_local_region_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    return lines

def parse_local_region(lines):
    row_start_y = 0
    row_height = 0
    rows = []
    cells = []

    for line in lines:
        parts = line.strip().split()
        if not parts:
            continue

        if parts[0] == 'RowStartY':
            row_start_y = int(parts[1])
        elif parts[0] == 'RowHeight':
            row_height = int(parts[1])
        elif parts[0] == 'Row':
            row_id = int(parts[1])
            start_x = int(parts[2])
            end_x = int(parts[3])
            rows.append((row_id, start_x, end_x))
        else:
            cell_name = parts[0]
            lower_left_x = int(parts[1])
            lower_left_y = int(parts[2])
            width = int(parts[3])
            height = int(parts[4])
            cells.append((cell_name, lower_left_x, lower_left_y, width, height))

    return row_start_y, row_height, rows, cells

def draw_local_region(row_start_y, row_height, rows, cells, output_path):
    fig, ax = plt.subplots()

    # Track the min and max coordinates for setting the axes limits
    min_x = float('inf')
    max_x = float('-inf')
    min_y = float('inf')
    max_y = float('-inf')

    # Draw rows
    for row_id, start_x, end_x in rows:
        row_y = row_start_y + row_height * row_id
        rect = patches.Rectangle((start_x, row_y), end_x - start_x, row_height,
                                 linewidth=0.1, edgecolor='green', facecolor='none', alpha=0.5)
        ax.add_patch(rect)

        min_x = min(min_x, start_x)
        max_x = max(max_x, end_x)
        min_y = min(min_y, row_y)
        max_y = max(max_y, row_y + row_height)

    # Draw cells
    for cell_name, lower_left_x, lower_left_y, width, height in cells:
        rect = patches.Rectangle((lower_left_x, lower_left_y), width, height,
                                 linewidth=0.1, edgecolor='blue', facecolor='none', alpha=0.5)
        ax.add_patch(rect)

        min_x = min(min_x, lower_left_x)
        max_x = max(max_x, lower_left_x + width)
        min_y = min(min_y, lower_left_y)
        max_y = max(max_y, lower_left_y + height)

    # Set axes limits
    print(min_x, max_x, min_y, max_y)
    ax.set_xlim(min_x - 10, max_x + 10)
    ax.set_ylim(min_y - 10, max_y + 10)

    ax.set_aspect('equal')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    plt.title('Local Region Visualization')

    # Save the plot as an image with higher DPI
    plt.savefig(output_path, dpi=1200)
    plt.close()

def main():
    if len(sys.argv) != 3:
        print("Usage: python draw_local_region.py <local_region.txt> <output_image.png>")
        sys.exit(1)

    file_path = sys.argv[1]
    output_path = sys.argv[2]
    lines = read_local_region_file(file_path)
    row_start_y, row_height, rows, cells = parse_local_region(lines)
    draw_local_region(row_start_y, row_height, rows, cells, output_path)

if __name__ == "__main__":
    main()
