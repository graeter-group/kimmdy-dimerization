ff_file = open("TdT_React_FF.top")
kimmdy_file = open("TdT_React_FF.top")
keyword = "pairs"
toggle_dih = False
dihedrals_ff = []
for line in ff_file:
    if toggle_dih and line.find("[") != -1:
        break
    if toggle_dih:
        rs = [a for a in line.strip("\n").split(" ") if a != ""]
        dihedrals_ff.append(rs)
    if line.find(f"[ {keyword} ]") != -1 and not toggle_dih:
        toggle_dih = True
        continue
    else:
        continue

toggle_dih = False
dihedrals_kimmdy = []
for line in kimmdy_file:
    if toggle_dih and line.find("[") != -1:
        break
    if toggle_dih:
        rs = [a for a in line.strip("\n").split(" ") if a != ""]
        dihedrals_kimmdy.append(rs)
    if line.find(f"[ {keyword} ]") != -1 and not toggle_dih:
        toggle_dih = True
        continue
    else:
        continue

for dih in dihedrals_kimmdy:
    if dih not in dihedrals_ff:
        print(dih)

for dih in dihedrals_ff:
    if dih not in dihedrals_kimmdy:
        print(dih)