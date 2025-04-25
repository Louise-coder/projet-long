def read(pdb_file: str, pred_start, exp_start):
    res_num = exp_start
    old_res_num = str(pred_start)
    with open("scripts/corrected.pdb", "w") as f_corr:
        with open("scripts/to_correct.pdb", "r") as f_to_corr:
            for line in f_to_corr:
                if line[77:78] == "H":
                    continue
                if line[22:26] != old_res_num:
                    res_num += 1
                    old_res_num = line[22:26]
                corrected_line = (
                    line[:22] + str(res_num) + line[26:72] + "A" + line[73:]
                )
                f_corr.write(corrected_line)


if __name__ == "__main__":
    pred_start = 1213
    exp_start = 1235
    read("scripts/to_correct.pdb", pred_start, exp_start)
