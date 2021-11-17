def DynProg(E_annual, YearsNum, E_opt, wells, r, discr_E):

    max_edge = math.sqrt(E_opt[0]**2 + E_opt[0]**2)
    dE = max_edge / discr_E
    print("dE = ", dE)



    DP_points = []

    point_last = Point(E_opt, wells)
    DP_points.insert(0, [point_last])
    for i, w in enumerate(wells):
        print("E opt ", str(i + 1), " = ", DP_points[-1][0].E[i], "\tQ opt ", str(i + 1), " = ", DP_points[-1][0].Q[i])


    # go back
    for main_iteration_counter in range(YearsNum - 1, -1, -1):

        #########################
        # 1) create a new points

        # cumulative external injection
        E_current = sum(E_annual[:main_iteration_counter])

        # find point according to E1 + E2 = E_current
        E_edge = [[0.0, 0.0], [0.0, 0.0]]

        for i,w in enumerate(wells):
            if E_current > E_opt[i]:
                E_edge[i][i] = E_opt[i]
                E_edge[i][(i+1)%2] = E_current - E_opt[i]
            else:
                E_edge[i][i] = E_current
                E_edge[i][(i+1)%2] = 0.0

        #print(E_edge)

        delta_edge = math.sqrt((E_edge[0][0] - E_edge[1][0])**2 + (E_edge[0][1] - E_edge[1][1])**2)
        #print("de = ", delta_edge)
        current_points_number = round(delta_edge / dE + 0.5) + 1
        #print("cpn = " , current_points_number)
        if delta_edge != 0.0:
            dE_current = delta_edge / (current_points_number - 1)
        else:
            dE_current = 0.0


        # define Points for current iteration
        current_points = []
        dr = [E_edge[1][0] - E_edge[0][0], E_edge[1][1] - E_edge[0][1]]
        for p in range(current_points_number):
            if current_points_number != 1:
                current_points.append(Point([E_edge[0][0] + p * dr[0] / (current_points_number - 1), E_edge[0][1] + p * dr[1] / (current_points_number - 1)], wells) )
            else:
                current_points.append(Point([E_edge[0][0], E_edge[0][1]], wells))

        DP_points.insert(0, current_points)

        #########################
        # 2) compute optimal value
        for point in DP_points[0]:

            if point.optimal_value != -1:
                print('error in opt value\n')
                exit(1)
            if point.optimal_pointer != -1:
                print('error in opt ptr\n')
                exit(1)

            for j, next_point in enumerate(DP_points[1]):
                if ItCanBeTheNextPoint(point, next_point):
                    value = 1. / (1. + r)**(main_iteration_counter + 0.5) * sum( [ (next_point.Q[i] - point.Q[i]) for i in range(len(wells))] ) + next_point.optimal_value
                    if value > point.optimal_value:
                        point.optimal_value = value
                        point.optimal_pointer = j



    f, ax = plt.subplots(1)
    #ax.set_aspect(1)
    #ax.set_xlim(0, sum(E_annual))
    #ax.set_ylim(0, sum(E_annual))
    plt.xticks(np.arange(0., sum(E_annual), 0.1))
    plt.yticks(np.arange(0., sum(E_annual), 0.1))
    ax.set_xlim(0, 0.75)
    ax.set_ylim(0, 0.4)

    plt.plot([E_opt[0], E_opt[0]], [0, E_opt[1]], 'k--')
    plt.plot([0, E_opt[0]], [E_opt[1], E_opt[1]], 'k--')
    #print(len(DP_points))
    for iteration in DP_points:
        #print("\t", len(iteration))
        for point in iteration:
            #print(point.E)
            plt.plot([point.E[0]], [point.E[1]], color='k', marker='.')

    # plot optimal way
    # go straight
    current_pointer = DP_points[0][0].optimal_pointer
    r0 = [DP_points[0][0].E[0], DP_points[0][0].E[1]]
    # print(current_pointer)
    for main_iteration_counter in range(1, YearsNum+1):
        # print(DP_points[main_iteration_counter][current_pointer].optimal_pointer)
        r1 = [DP_points[main_iteration_counter][current_pointer].E[0],
                DP_points[main_iteration_counter][current_pointer].E[1]]
        plt.plot([r0[0], r1[0]], [r0[1], r1[1]], 'r-')
        #plt.arrow(r0[0], r0[1], r1[0] - r0[0], r1[1] - r0[1], color = 'r', ls='-', lw = 2)
        current_pointer = DP_points[main_iteration_counter][current_pointer].optimal_pointer
        r0 = r1

    #plt.savefig('figure3.pdf', format='pdf')
    #plt.show()
    plt.show()


################################################################################

def ItCanBeTheNextPoint(point, next_point):
    for i in range(len(point.E)):
        if next_point.E[i] < point.E[i]:
            return False

    return True
