import sys
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def main():
    # 檢查命令列參數
    if len(sys.argv) != 3:
        print("使用方法: python draw_lg.py <輸入檔案路徑> <輸出圖形路徑>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    # 讀取輸入檔案內容
    with open(input_file, 'r', encoding='utf-8') as f:
        data = f.read()

    # 解析資料
    lines = data.strip().split('\n')

    die_size = None
    cells = []
    placement_rows = []

    for line in lines:
        tokens = line.strip().split()
        if not tokens:
            continue
        if tokens[0] == 'Alpha' or tokens[0] == 'Beta':
            continue  # 忽略這些行
        elif tokens[0] == 'DieSize':
            die_size = list(map(int, tokens[1:]))
        elif tokens[0] == 'PlacementRows':
            start_x = int(tokens[1])
            start_y = int(tokens[2])
            site_width = int(tokens[3])
            site_height = int(tokens[4])
            total_sites = int(tokens[5])
            placement_rows.append({
                'x': start_x,
                'y': start_y,
                'width': site_width * total_sites,
                'height': site_height
            })
        else:
            # 解析 cell 資料
            cell_name = tokens[0]
            lower_left_x = int(tokens[1])
            lower_left_y = int(tokens[2])
            width = int(tokens[3])
            height = int(tokens[4])
            fix_status = tokens[5]
            cells.append({
                'name': cell_name,
                'x': lower_left_x,
                'y': lower_left_y,
                'width': width,
                'height': height,
                'status': fix_status
            })

    # 檢查是否有 DieSize 資料
    if not die_size:
        print("未找到 DieSize 資料，請確認輸入檔案格式正確。")
        sys.exit(1)

    # 設定圖像解析度
    dpi = 600
    linewidth_pt = 2*(72 / dpi)  # 對應於 1 像素的 linewidth

    # 開始繪圖
    fig, ax = plt.subplots(figsize=(10, 10))

    # 繪製 DieSize 邊框
    die_rect = patches.Rectangle(
        (die_size[0], die_size[1]),
        die_size[2] - die_size[0],
        die_size[3] - die_size[1],
        edgecolor='#11111b',
        facecolor='none',
        linewidth=2,
        alpha=1
    )
    ax.add_patch(die_rect)

    # 定義顏色對應
    color_map = {
        'FIX': 'red',
        'NOTFIX': 'blue'
    }

    # 繪製 PlacementRows，邊框為綠色
    for row in placement_rows:
        rect = patches.Rectangle(
            (row['x'], row['y']),
            row['width'],
            row['height'],
            linewidth=linewidth_pt/4,
            edgecolor='#1e1e2e',
            facecolor='orange',
            alpha=0.3
        )
        ax.add_patch(rect)

    # 繪製 cells，僅繪製邊框，邊框粗細為 1 像素
    for cell in cells:
        rect = patches.Rectangle(
            (cell['x'], cell['y']),
            cell['width'],
            cell['height'],
            linewidth=linewidth_pt/4,
            edgecolor='#1e1e2e',
            facecolor=color_map.get(cell['status'], '#4c4f69'),
            alpha=0.5
        )
        ax.add_patch(rect)

    # 設定圖形範圍和標籤
    ax.set_xlim(die_size[0], die_size[2])
    ax.set_ylim(die_size[1], die_size[3])
    ax.ticklabel_format(style='plain')  # 禁用科學記號
    ax.set_xticks(range(0, die_size[0] + 1, max(1, die_size[0] // 10)))  # 使用整數標籤
    ax.set_yticks(range(0, die_size[1] + 1, max(1, die_size[1] // 10)))  # 使用整數標籤
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel('X Axis')
    plt.ylabel('Y Axis')
    plt.title('Cell and PlacementRows Visualization')
    # plt.legend(handles=[
    #     patches.Patch(edgecolor='blue', facecolor='none', label='NOTFIX'),
    #     patches.Patch(edgecolor='red', facecolor='none', label='FIX'),
    #     patches.Patch(edgecolor='green', facecolor='none', label='PlacementRows')
    # ])

    # 儲存圖形到指定的輸出路徑，確保圖像解析度
    # Show x=821370~1231770, y=991200~993600
    plt.savefig(output_file, dpi=dpi)
    plt.close()
    print(f"圖形已儲存到 {output_file}")

if __name__ == "__main__":
    main()
