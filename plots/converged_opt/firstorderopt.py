import niceplots
import matplotlib.pyplot as plt

data = []
iterations = []
optimality = []
feasibility = []
obj_function = []
with open("firstorderopt.txt", "r") as f:
    for line in f:
        # Split line by whitespace
        line_data = line.split()
        # Extract major, minors, step, nCon, feasible, optimal, merit function, and penalty columns
        major = int(line_data[0])
        minors = int(line_data[1])
        step = (line_data[2])
        nCon = (line_data[3])
        feasible = (line_data[4])
        feasible = feasible.strip("(")
        feasible = feasible.strip(")")
        feasible1, feasible2 = feasible.split("E")
        feasible = float(feasible1) * 10**(int(feasible2))
        optimal = (line_data[5])
        optimal = optimal.strip("(")
        optimal = optimal.strip(")")
        optimal1, optimal2 = optimal.split("E")
        optimal = float(optimal1) * 10**(int(optimal2))
        merit_function = (line_data[6])
        merit_function = merit_function.strip("(")
        merit_function = merit_function.strip(")")
        merit_function1, merit_function2 = merit_function.split("E")
        merit_function = float(merit_function1) * 10**(int(merit_function2)-1)
        # Append data as a tuple to the list
        data.append((major, minors, step, nCon, feasible, optimal, merit_function))
        iterations.append(major)
        optimality.append(optimal)
        feasibility.append(feasible)
        obj_function.append(merit_function)

niceplots.setRCParams()
_, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(8, 10))
ax1.plot(iterations[:40], optimality[:40])
ax1.set_ylabel('Optimality')
ax1.set_ylim((-1e-2, 5e-2))
ax2.plot(iterations[:40], feasibility[:40])
ax2.set_ylabel('Feasibility')
ax2.set_ylim((-1e-2, 6e-2))
ax3.plot(iterations[:40], obj_function[:40])
ax3.set_ylabel('Objective function')
ax3.set_xlabel('Iterations')
plt.savefig('optimisation_040.png')
plt.savefig('optimisation_040.eps', format='eps')

_, (ax4,ax5,ax6) = plt.subplots(3,1,figsize=(8, 10))
ax4.plot(iterations[40:], optimality[40:])
ax4.set_ylabel('Optimality')
ax4.set_ylim((-1e-6, 5e-4))
ax5.plot(iterations[40:], feasibility[40:])
ax5.set_ylabel('Feasibility')
ax5.set_ylim((-1e-8, 6e-7))
ax6.plot(iterations[40:], obj_function[40:])
ax6.set_ylabel('Objective function')
ax6.set_xlabel('Iterations')
niceplots.adjust_spines(ax1, outward=True)
niceplots.adjust_spines(ax2, outward=True)
niceplots.adjust_spines(ax3, outward=True)
niceplots.adjust_spines(ax4, outward=True)
niceplots.adjust_spines(ax5, outward=True)
niceplots.adjust_spines(ax6, outward=True)
plt.savefig('optimisation_4080.png')
plt.savefig('optimisation_4080.eps', format='eps')