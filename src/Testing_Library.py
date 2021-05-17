from Generate_Seq import SimulateSeqCD

domain_choices_var = [2, 4, 6, 8, 10, 20, 30, 40, 50]
domain_choices_eq = [1, 2, 4, 5, 10, 20, 25, 40, 50, 80, 100]
sim = SimulateSeqCD("C:/Users/RackS/Documents/Isochores100/")
trials = 100

#generate equal-length test sequences:
for d in domain_choices_eq:
    for i in range(0, trials):
        name = "E" + (str(d).zfill(3)) + "_" + (str(i+1).zfill(3))
        sim.generate_equal(num_domains=d, id=name)
    print()

#generate variable-length test sequences:
for d in domain_choices_var:
    for i in range(0, trials):
        name = "V" + (str(d).zfill(3)) + "_" + (str(i+1).zfill(3))
        sim.generate_variable(num_domains=d, id=name)
    print()