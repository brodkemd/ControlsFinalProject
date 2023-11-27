import numpy as np
import logging, sys, os, json
sys.path.append(".")
from tools.rotations import Euler2Quaternion
from tools.toString import arrToStr
from itertools import permutations 
import parameters.simulation_parameters as SIM
from viewers.dataPlotter                import DataPlotter
from dynamics.dynamics                  import Dynamics
from dynamics.forcesMoments             import ForcesMoments
from controller.fullStateFeedback       import FullStateFeedBack



def descent():
    cwd = os.path.dirname(__file__)
    # Configure the logging system
    logging.basicConfig(level=logging.DEBUG,
                        format='descent, %(levelname)s: %(message)s')

    # Create a logger
    logger = logging.getLogger(__name__)

    with open(os.path.join(cwd, 'descent.log'), 'w') as f:
        pass
    # Create a handler for writing log messages to a file
    file_handler = logging.FileHandler(os.path.join(cwd, 'descent.log'))
    file_handler.setLevel(logging.DEBUG)  # Set the level for the file handler

    # Create a handler for displaying log messages on the console
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)  # Set the level for the console handler

    # Create a formatter
    formatter = logging.Formatter('%(levelname)s: %(message)s')

    # Set the formatter for both handlers
    file_handler.setFormatter(formatter)
    console_handler.setFormatter(formatter)

    # Add the handlers to the logger
    logger.addHandler(file_handler)
    logger.addHandler(console_handler)

    plotter = DataPlotter(width = SIM.fig_width, height = SIM.fig_height, plot = False, interactive = False)

    phi   = 0
    theta = np.deg2rad(0) # np.deg2rad(90)
    psi   = 0
    state = np.array([
        -10, # p_n
         0, # p_e
        -11000, # p_d
         0, # u
         0, # v
         182.45, # w, M=0.8 for entry
         Euler2Quaternion(phi, theta, psi).item(0), # e_0
         Euler2Quaternion(phi, theta, psi).item(1), # e_1
         Euler2Quaternion(phi, theta, psi).item(2), # e_2
         Euler2Quaternion(phi, theta, psi).item(3), # e_3
         0, # p
         0, # q
         0 # r
    ], dtype=float)
    perm = list(permutations(5*[10]+5*[1]+5*[0.1]+5*[0.01], 5))

    data = {}
    perm_count = 0
    for i in range(len(perm)):
        print(i, "/", len(perm), end="    \n")
        data[perm_count] = {}
        try:
            controller = FullStateFeedBack(compute_gains=True)
            controller.descentCP.omega_ns = np.array(perm[i], dtype=float)
            data[perm_count]["omega_ns"] = list(perm[i])
            logger.info('Running with omega_ns:' + arrToStr(controller.descentCP.omega_ns, num_format=":.3f"))

            controller.descentCP.computePoles()
            controller.descentCP.generateGains(compute_gains=True)
            dynamics   = Dynamics(state.copy())
            forces     = ForcesMoments()
        
        except KeyboardInterrupt:
            # catches user interrupts, stops wordy error messages that occur due a usual user interruption
            print("\nCaught keyboard interrupt, Exiting ...")
            break

        except Exception as e:
            logger.error("  "+str(e))
            data[perm_count]["succeed"] = False
            perm_count+=1
            continue

        try:
            print("\nProgress:")
            t = SIM.start_time
            count = 0
            while t < SIM.end_time:
                # setting next time to make a plot, exits inner loop at that time
                t_next_plot = t + SIM.ts_plot

                # loop that runs the dynamics
                while t < t_next_plot:
                    F_E, F_cp, tau, x_r = controller.update(state)
                    u                   = forces.update(state, F_E, F_cp, tau)
                    y, crash, landed    = dynamics.update(u)
                    t += SIM.ts_simulation
                    if crash or landed: break

                # plotting the state variables and forces with respect to time
                response = [[y.item(i), x_r.item(i)] for i in range(len(y))]
                plotter.update(t, response, F_E.tolist() + F_cp.tolist() + tau.tolist())

                # -------increment time-------------
                count+=1

                if crash:
                    print("Detected crash :(\nexiting simulation loop")
                    break
                if landed:
                    print("Detected landing, yay!!!!!\nexiting simulation loop")
                    break

                # progresses bar with stats
                #print(round(t, 6), "/", SIM.end_time, round(np.rad2deg(Quaternion2Euler(y[6:10])[1]),2), end="      \r")
                print(round(t, 6), "/", SIM.end_time, end="      \r")
            print()
            logger.info("  Completed successfully")
            data[perm_count]["succeed"] = True

        except KeyboardInterrupt:
            # catches user interrupts, stops wordy error messages that occur due a usual user interruption
            print("\nCaught keyboard interrupt, Exiting ...")
            break

        except Exception as e:
            # saves data if there is an error
            # traceback.print_exc()
            logger.error("  "+str(e))
            data[perm_count]["succeed"] = False
            logger.info("  Dumping collected data")
        plotter.saveData(os.path.join(cwd, f"{perm_count}_data.json"))
        perm_count+=1
    
    with open(os.path.join(cwd,"descent_info.json"), 'w') as f:
        f.write(json.dumps(data, indent=2))

descent()