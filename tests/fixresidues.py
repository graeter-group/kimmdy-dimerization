
gro_file_name = "DoubleStrand_NPT.gro"
top_file_name = "DoubleStrand_NPT.top"
gro_file_name_fixed = "DoubleStrand_NPT_fixed.gro"
top_file_name_fixed = "DoubleStrand_NPT_fixed.top"

gro_file = open(gro_file_name)
gro_file_out = open(gro_file_name_fixed, "w")
cur_res = 1
last_res = 1
found_idx = 0

for line in gro_file:
    if line.find("D") == -1:
        gro_file_out.write(line)
    else:
        line_par = line.strip("\n").split(" ")
        if line_par[3].find("D") == -1:
            res_entry = line_par[4]
            found_idx = 4
        else:
            res_entry = line_par[3]
            found_idx = 3
        residue = int(res_entry.split("D")[0])
        if last_res != residue:
            if residue < cur_res:
                cur_res += 1
            else:
                cur_res = residue

        line_par[found_idx] = (str(cur_res) + "D" + res_entry.split("D")[1])
        if found_idx == 4 and cur_res > 9:
            gro_file_out.write(" ".join(line_par[1:]) + "\n")
        else:
            gro_file_out.write(" ".join(line_par) + "\n")

        last_res = residue

gro_file_out.close()
gro_file.close()

top_file = open(top_file_name)
top_file_out = open(top_file_name_fixed, "w")

cur_res = 1
last_res = 1
found_idx = 0

for line in top_file:
    if line.find("D") == -1 or (line.find(";") != -1  and line.find("qtot") == -1) or line.find("DNA") != -1:
        top_file_out.write(line)
    else:
        line_par = line.strip("\n").split(" ")
        line_par_short = [a for a in line_par if a != ""]
        residue = int(line_par_short[2])
        if last_res != residue:
            if residue < cur_res:
                cur_res += 1
            else:
                cur_res = residue

        for i in range(14,23):
            if line_par[i] != "":
                try:
                    res = int(line_par[i])
                    found_idx = i
                except:
                    pass

        line_par[found_idx] = str(cur_res)
        top_file_out.write(" ".join(line_par) + "\n")
        last_res = residue


top_file_out.close()
top_file.close()