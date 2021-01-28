import numpy as np
import matplotlib.pyplot as plt
from random import random
import random as rand
import math

# Simulation of interval ALOHA algorithm

# MODE:
# a - The user randomly selects the window number for transmission from the interval {1 2, ... , W}.
# b - The first transmission of a message occurs at the beginning of the next window, during 
#     subsequent transmissions of the message, the user accidentally
#     selects a window number from an interval {1 2, ... , W}.
# s - comparison of simulation results for M = 1 with the synchronous system M|D|1.
MODE = 'a'

class Message:
    slots_left_before_transfer = -1
    time_of_appearance_in_system = -1
    def __init__(self, time_of_appearance_in_system):
        self.time_of_appearance_in_system = time_of_appearance_in_system
    def tick(self):
        self.slots_left_before_transfer = self.slots_left_before_transfer - 1

def generate_stream(lambd, num_windows, num_of_users):
    when_user_message_appears_in_system = [] 
    # Poisson Distribution for M users
    # for i in range(num_of_users):
    #     T = [] 
    #     while sum(T) < num_windows: # As long as the sum T is not more than the total simulation time,
    #         t = ((-num_of_users / lambd) * (math.log(random()))) # generate t exponentially 
    #         T.append(t) # and add to the general timeline T
    #     T.pop() # Removing the last t, since it is exactly greater than T
    #     T = np.cumsum(T) # Change T to temporary points (instead of segments)
    #     when_user_message_appears_in_system.append(T) 

    for _ in range(num_of_users):
        tau = [((-num_of_users/lambd)*math.log(random())) for _ in range(num_windows)]
        t = np.cumsum(tau)
        when_user_message_appears_in_system.append(t)
    
    return when_user_message_appears_in_system

def synchronous_modeling (stream1, num_windows):
    stream = []
    num_requests_per_windows = [0]*num_windows 

    for i in range(len(stream1)):
        for j in range(len(stream1[i])):
            stream.append(stream1[i][j])
    stream.sort()

    for j in range(len(stream)):
        window = math.trunc(stream[j])
        if  window < num_windows:
            num_requests_per_windows[window] += 1

    num_requests_now = 0 
    num_requests_quit = 0
    current_request = 0
    D = 0
    N = 0

    for i in range(num_windows): 
        if num_requests_now != 0: 
            num_requests_now -= 1 # send request
            D += i + 1 - stream[current_request]
            current_request += 1 
            num_requests_quit += 1
        num_requests_now += num_requests_per_windows[i] # add requests that appeared in this window 
        N += num_requests_now

    if num_requests_quit == 0:
        M_D = 0
    else:
        M_D = D / num_requests_quit
    M_N = N / num_windows 
    L_out = num_requests_quit / num_windows

    return M_D, M_N, L_out

def interval_aloha(stream, num_windows, window_size_W):
    
    num_of_users = len(stream)
    D = 0
    N = 0
    num_requests_quit = 0
    user_message_queue = [[] for _ in range(num_of_users)]
    count_collision = 0

    for current_window in range(num_windows):
        transmitted_nodes = []

        # polling all users in current window
        for current_user in range(num_of_users): 
            for elem in stream[current_user]: # how many messages does user want to send in current window          
                if math.trunc(elem) == current_window:
                    user_message_queue[current_user].append(Message(elem))
                    stream[current_user] = np.delete(stream[current_user], 0, 0) # delete first elem of array stream[i]
                elif math.trunc(elem) > current_window:
                    break

            # take first element of queue
            if len(user_message_queue[current_user]) == 1 and MODE == 'b':
                user_message_queue[current_user][0].slots_left_before_transfer = 1

            if len(user_message_queue[current_user]) != 0:
                if user_message_queue[current_user][0].slots_left_before_transfer < 0:
                    user_message_queue[current_user][0].slots_left_before_transfer = rand.randint(1, window_size_W)
                else:
                    user_message_queue[current_user][0].tick()
                    
                if user_message_queue[current_user][0].slots_left_before_transfer == 0:
                    transmitted_nodes.append(current_user)
        
        # success
        if (len(transmitted_nodes) == 1): 
            num_requests_quit += 1 
            D += (current_window + 1 - user_message_queue[transmitted_nodes[0]][0].time_of_appearance_in_system)
            # throw coin to pick a window
            if len(user_message_queue[transmitted_nodes[0]]) > 1 and user_message_queue[transmitted_nodes[0]][1].slots_left_before_transfer < 0:
                    user_message_queue[transmitted_nodes[0]][1].slots_left_before_transfer = rand.randint(1, window_size_W)
            # remove a message from the queue
            user_message_queue[transmitted_nodes[0]].pop(0)

        # collision
        if (len(transmitted_nodes) > 1): 
            count_collision += 1
            for j in transmitted_nodes:
                user_message_queue[j][0].slots_left_before_transfer = rand.randint(1, window_size_W)
        
        for i in range(num_of_users):
                if len(user_message_queue[i]) != 0:
                    N += len(user_message_queue[i])
        
    M_D = D / num_requests_quit
    M_N = N / num_windows 
    # print("num of collision: ", count_collision)
    L_out = num_requests_quit/num_windows
    return M_D, M_N, L_out



def main():
    num_windows = 10000
    M = 5 # number of users
    W = 9 # transmission interval
    lambdas = np.arange(0.1, 1.0, 0.05)
    list_D = []
    list_N = []
    lambdas_out_list = []
    list_D_s = []
    list_N_s = []
    lambdas_out_list_s = []
    list_D_s_theor = []
    list_N_s_theor = []

    if MODE == 's':
        M = 1
        W = 1

    P = 2/(W+1)
    lambd_critical = M*P*((1-P)**(M-1))

    for l in lambdas:
        print("λ = ", round(l,2))
        stream = generate_stream(l, num_windows, M)
        if MODE == 's':
            # Synchronous system
            M_D_s, M_N_s, lambda_out_s = synchronous_modeling(stream, num_windows)
            lambdas_out_list_s.append(lambda_out_s)
            list_D_s.append(M_D_s)
            list_N_s.append(M_N_s)
            print("modeling M[D] (Synchr) = ", M_D_s)
            print("modeling M[N] (Synchr) = ", M_N_s)
        # Interval aloha
        M_D, M_N, lambda_out  = interval_aloha(stream, num_windows, W)
        print("modeling M[D] (Interval) = ", M_D)
        print("modeling M[N] (Interval) = ", M_N)
        lambdas_out_list.append(lambda_out)
        list_D.append(M_D)
        list_N.append(M_N)
        
        # Synchronous system (Theor)
        # list_D_s_theor.append((3 - 2*l) / (2*(1 - l)))
        # list_N_s_theor.append((l*(2 - l)) / (2*(1 - l)))
        # print("modeling M[N] (Synchr theor) = ", (l*(2 - l)) / (2*(1 - l)))


    plt.figure(1)
    plt.xlabel('lambda')
    plt.axvline(lambd_critical, label='lambda critical = ' + str(round(lambd_critical,2)), color='r')
    plt.plot(lambdas, list_D, label ='M_D (Interval)')
    plt.plot(lambdas, list_N, label ='M_N (Interval)')
    if MODE == 's':
        plt.plot(lambdas, list_D_s,label ='M_D (Synchr)', linestyle = '--')
        plt.plot(lambdas, list_N_s,label ='M_N (Synchr)', linestyle = '--')
        # plt.plot(lambdas, list_D_s_theor,label='M_D (Synchr theor)')
        # plt.plot(lambdas, list_N_s_theor,label='M_N (Synchr theor)')
    plt.legend()
    plt.show()
    plt.figure(2)
    plt.xlabel('lambda in')
    plt.ylabel('lambda out')
    plt.axhline(lambd_critical, label='lambda critical = ' + str(round(lambd_critical,2)), color='r')
    plt.plot(lambdas, lambdas_out_list, label='lambda out (Interval)', color = 'darkmagenta')
    if MODE == 's':
        plt.plot(lambdas, lambdas_out_list_s, label='lambda out (Synchr)', linestyle = '--', color='y')
    plt.legend()
    plt.show()
    plt.close()
    
if __name__ == "__main__":
	main()