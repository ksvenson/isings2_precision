all: gen_grp_o3 gen_grp_o4

.phony: gen_grp_o3 gen_grp_o4 clean

clean:
	$(RM) -rf elem* lattice* site_g*

gen_grp_o3: gen_grp_o3.py grp_o3.py quat.py
	python3 gen_grp_o3.py

gen_grp_o4: gen_grp_o4.py grp_o4.py quat.py
	python3 gen_grp_o4.py
