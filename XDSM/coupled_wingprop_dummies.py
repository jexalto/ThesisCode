from pyxdsm.XDSM import (
    XDSM,
    OPT,
    SUBOPT,
    SOLVER,
    DOE,
    IFUNC,
    FUNC,
    GROUP,
    IGROUP,
    METAMODEL,
    LEFT,
    RIGHT,
)

# Change `use_sfmath` to False to use computer modern
x = XDSM(use_sfmath=True)
optlst = [r"l^*, \mathbf{\phi}^*, c_r^*, c_t^*", r"B^*, R^*, \mathbf{\theta}^*", r"\mathbf{c}^*, J^*, c^*"]
optlst0 = [r"l^0, \mathbf{\phi}^0, c_r^0, c_t^0", r"B^0, R^0, \mathbf{\theta}^0", r"\mathbf{c}^0, J^0, c^0"]
con1lst = [r"\text{Thrust}", r"\text{constraint}"]
con2lst = [r"\text{Lift}", r"\text{constraint}"]

x.add_system("opt", OPT, r"\text{Optimizer (SNOPT)}")
x.add_system("MDA", SOLVER, r"\text{MDA coordinator}")
x.add_system("propeller", FUNC, r"\text{Propeller}")
x.add_system("slipstream", FUNC, r"\text{Slipstream model}")
x.add_system("wing", FUNC, r"\text{Wing}")
x.add_system("performance", FUNC, r"\text{Performance}")
x.add_system("F", IFUNC, r"\text{Objective function}")
x.add_system("G", IFUNC, con1lst)
x.add_system("G2", IFUNC, con2lst)

lst = [r"B, R, \mathbf{\theta}", r"\mathbf{c}, J, c"]
x.connect("opt", "propeller", lst)
x.connect("opt", "wing", r"l, \mathbf{\phi}, c_r, c_t")
x.connect("MDA", "wing", r'\hat{W}_{fuel}')
x.connect("MDA", "propeller", r'\hat{\alpha}')
x.connect("propeller", "slipstream", r"\mathbf{v_j}")
x.connect("propeller", "performance", r"Power")
x.connect("slipstream", "wing", r"\mathbf{G}, mesh")
x.connect("wing", "performance", r"Lift, Drag, W_{str}")
x.connect("wing", "MDA", r"\alpha")
x.connect("performance", "MDA", r"W_{fuel}")

x.connect("propeller", "F", r"Power")
x.connect("propeller", "G", r"Thrust")
x.connect("wing", "G", r"Drag")
x.connect("wing", "G2", r"Lift, W_{structural}")
x.connect("performance", "G2", r"W_{fuel}")

x.connect("F", "opt", "Power")
x.connect("G", "opt", r"Thrust=Drag")
x.connect("G2", "opt", r"Lift=Weight")

x.add_process(["opt", "MDA", "propeller", "slipstream", "wing", "performance", "F", "opt"],
    arrow=True,)
x.add_process(["propeller", "performance"],
    arrow=True,)
x.add_process(["propeller", "F"],
    arrow=True,)
x.add_process(["propeller", "G"],
    arrow=True,)
x.add_process(["wing", "MDA"],
    arrow=True,)
x.add_process(["wing", "G", "opt"],
    arrow=True,)
x.add_process(["wing", "G2"],
    arrow=True,)
x.add_process(["performance", "G2", "opt"],
    arrow=True,)
x.add_process(["performance", "MDA"],
    arrow=True,)
x.add_process(["MDA", "wing"],
    arrow=True,)

x.add_input("MDA", "W_f^0")
x.add_input("opt", optlst0)

# x.add_output("opt", optlst, side=LEFT)
# x.add_output("propeller", r"\mathbf{v_j}^*, P_{in}^*, T^*", side=LEFT)
# x.add_output("slipstream", r"\mathbf{G}^*", side=LEFT)
# x.add_output("wing", r"C_L^*, C_D^*, W_{str}^*", side=LEFT)
# x.add_output("performance", r"W_f^*", side=LEFT)
# x.add_output("F", r"f^*", side=LEFT)
# x.add_output("G", r"g_1^*", side=LEFT)
# x.add_output("G2", r"g_2^*", side=LEFT)

x.write("coupled_wingprop_fordummies")
