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
x.add_system("helix", FUNC, r"\text{HELIX}")
x.add_system("slipstream", FUNC, r"\text{Slipstream model}")
x.add_system("OAS", FUNC, r"\text{EnhancedOpenAeroStruct}")
x.add_system("performance", FUNC, r"\text{Performance}")
x.add_system("F", IFUNC, r"\text{Objective function}")
x.add_system("G", IFUNC, con1lst)
x.add_system("G2", IFUNC, con2lst)

lst = [r"B, R, \mathbf{\theta}", r"\mathbf{c}, J, c"]
x.connect("opt", "helix", lst)
x.connect("opt", "OAS", r"l, \mathbf{\phi}, c_r, c_t")
x.connect("MDA", "OAS", r'\hat{W}_f')
x.connect("MDA", "helix", r'\hat{\alpha}')
x.connect("helix", "slipstream", r"\mathbf{v_j}")
x.connect("helix", "performance", r"P_{in}")
x.connect("slipstream", "OAS", r"\mathbf{G}, mesh")
x.connect("OAS", "performance", r"C_L, C_D, W_{str}")
x.connect("OAS", "MDA", r"\alpha")
x.connect("performance", "MDA", r"W_f")

x.connect("helix", "F", r"P_{in}")
x.connect("helix", "G", r"T")
x.connect("OAS", "G", r"C_D")
x.connect("OAS", "G2", r"C_L, W_{str}")
x.connect("performance", "G2", r"W_f")

x.connect("F", "opt", "f")
x.connect("G", "opt", r"g_1")
x.connect("G2", "opt", r"g_2")

x.add_process(["opt", "MDA", "helix", "slipstream", "OAS", "performance", "F", "opt"],
    arrow=True,)
x.add_process(["helix", "performance"],
    arrow=True,)
x.add_process(["helix", "F"],
    arrow=True,)
x.add_process(["helix", "G"],
    arrow=True,)
x.add_process(["OAS", "MDA"],
    arrow=True,)
x.add_process(["OAS", "G", "opt"],
    arrow=True,)
x.add_process(["OAS", "G2"],
    arrow=True,)
x.add_process(["performance", "G2", "opt"],
    arrow=True,)
x.add_process(["performance", "MDA"],
    arrow=True,)
x.add_process(["MDA", "OAS"],
    arrow=True,)

x.add_input("MDA", "W_f^0")
x.add_input("opt", optlst0)

x.add_output("opt", optlst, side=LEFT)
x.add_output("helix", r"\mathbf{v_j}^*, P_{in}^*, T^*", side=LEFT)
x.add_output("slipstream", r"\mathbf{G}^*", side=LEFT)
x.add_output("OAS", r"C_L^*, C_D^*, W_{str}^*", side=LEFT)
x.add_output("performance", r"W_f^*", side=LEFT)
x.add_output("F", r"f^*", side=LEFT)
x.add_output("G", r"g_1^*", side=LEFT)
x.add_output("G2", r"g_2^*", side=LEFT)

x.write("coupled_wingprop")
