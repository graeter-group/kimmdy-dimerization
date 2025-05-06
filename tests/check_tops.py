keywords = ["pairs", "angles", "bonds", "dihedrals", "atoms"]
morse_indices = [0, 8, 2, 3, 4, 5, 9, 10]
no_relax_indices = [0, 1, 2, 3, 4, 5, 6, 7]
indices = no_relax_indices

for keyword in keywords:
    ff_file = open("TdT_from_pdb2gmx.top")
    kimmdy_file = open("TdT_norelax.top")
    toggle_dih = False
    dihedrals_ff = []
    for line in ff_file:
        if line.find("vanishing") != -1:
            continue
        line = line.split(";")[0]
        if toggle_dih and line.find("[") != -1:
            break
        if toggle_dih:
            if keyword == "pairs":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:3]
            elif keyword == "dihedrals":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:5]
            elif keyword == "angles":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:4]
            elif keyword == "bonds":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:2]
            elif keyword == "atoms":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:8]
                rs[-2] = rs[-2][0:4]
            else:
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
        if line.find("vanishing") != -1:
            continue
        line = line.split(";")[0]
        if toggle_dih and line.find("[") != -1:
            break
        if toggle_dih:
            if keyword == "pairs":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:3]
            elif keyword == "dihedrals":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:5]
            elif keyword == "angles":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:4]
            elif keyword == "bonds":
                rs = [a for a in line.strip("\n").split(" ") if a != ""][0:2]
            elif keyword == "atoms":
                rs = [a for a in line.strip("\n").split(" ") if a != ""]
                rs = [rs[i] for i in indices]
                rs[-2] = rs[-2][0:4]
            else:
                rs = [a for a in line.strip("\n").split(" ") if a != ""]
            dihedrals_kimmdy.append(rs)
        if line.find(f"[ {keyword} ]") != -1 and not toggle_dih:
            toggle_dih = True
            continue
        else:
            continue

    for dih in dihedrals_kimmdy:
        if dih not in dihedrals_ff:
            print(keyword, dih)

    for dih in dihedrals_ff:
        if dih not in dihedrals_kimmdy:
            print(keyword, dih)

    ff_file.close()
    kimmdy_file.close()
