import numpy as np
import matplotlib.pyplot as plt
from random import random
import random as rand
import math

# Simulation of multichannel interval ALOHA algorithm
# Channel selection when slots_left_before_transfer (counter) = 0

class Message:
    slots_left_before_transfer = -1
    channel_number = -1
    time_of_appearance_in_system = -1
    def __init__(self, time_of_appearance_in_system):
        self.time_of_appearance_in_system = time_of_appearance_in_system
    def tick(self):
        self.slots_left_before_transfer = self.slots_left_before_transfer - 1

def generate_stream(lambd, num_windows, num_of_users):
    when_user_message_appears_in_system = [] 
    # Poisson Distribution for M users
    for _ in range(num_of_users):
        tau = [((-num_of_users/lambd)*math.log(random())) for _ in range(num_windows)]
        t = np.cumsum(tau)
        when_user_message_appears_in_system.append(t)

    return when_user_message_appears_in_system

def interval_aloha(stream, num_of_channels, num_windows, window_size_W):
    
    num_of_users = len(stream)
    D = 0
    N = 0
    num_requests_quit = 0
    user_message_queue = [[] for _ in range(num_of_users)]

    for current_window in range(num_windows):
        transmitted_nodes = [[] for _ in range(num_of_channels)]

        # polling all users in current window
        for current_user in range(num_of_users): 
            for elem in stream[current_user]: # how many messages does user want to send in current window          
                if math.trunc(elem) == current_window:
                    user_message_queue[current_user].append(Message(elem))
                    stream[current_user] = np.delete(stream[current_user], 0, 0) # delete first elem of array stream[i]
                elif math.trunc(elem) > current_window:
                    break

            if len(user_message_queue[current_user]) != 0:
                if user_message_queue[current_user][0].slots_left_before_transfer < 0:
                    user_message_queue[current_user][0].slots_left_before_transfer = rand.randint(1, window_size_W)
                else:
                    user_message_queue[current_user][0].tick()
                    
                if user_message_queue[current_user][0].slots_left_before_transfer == 0:
                    transmission_window = rand.randint(0, num_of_channels - 1)
                    user_message_queue[current_user][0].channel_number = transmission_window
                    transmitted_nodes[transmission_window].append(current_user)
        
        for current_channel in range(num_of_channels): 
            # success
            if (len(transmitted_nodes[current_channel]) == 1): 
                num_requests_quit += 1 
                D += (current_window + 1 - user_message_queue[transmitted_nodes[current_channel][0]][0].time_of_appearance_in_system)
                # throw coin to pick a window
                if len(user_message_queue[transmitted_nodes[current_channel][0]]) > 1 and user_message_queue[transmitted_nodes[current_channel][0]][1].slots_left_before_transfer < 0:
                        user_message_queue[transmitted_nodes[current_channel][0]][1].slots_left_before_transfer = rand.randint(1, window_size_W)
                # remove a message from the queue
                user_message_queue[transmitted_nodes[current_channel][0]].pop(0)

            # collision
            if (len(transmitted_nodes[current_channel]) > 1): 
                for j in transmitted_nodes[current_channel]:
                    user_message_queue[j][0].slots_left_before_transfer = rand.randint(1, window_size_W)
            
        for i in range(num_of_users):
                if len(user_message_queue[i]) != 0:
                    N += len(user_message_queue[i])
        
    M_D = D / num_requests_quit
    M_N = N / num_windows 
    L_out = num_requests_quit/num_windows
    return M_D, M_N, L_out


def main():
    num_windows = 10000
    M = 5 # number of users
    W = 3 # transmission interval
    K = 3 # number of channels
    
    P = 2/(W+1)
    lambd_critical = M*P*((1-P*(1/K))**(M-1))
    #lambd_critical = K*M*P*((1-P)**((M-1)/(K-1)))
    print('lambda critical = ' + str(round(lambd_critical,3)))

    lambdas = np.arange(0.1, 1.6, 0.05)
    list_D = []
    list_N = []
    lambdas_out_list = []

    for l in lambdas:
        print("λ = ", round(l,2))
        stream = generate_stream(l, num_windows, M)
        # Interval aloha
        M_D, M_N, lambda_out  = interval_aloha(stream, K, num_windows, W)
        print("modeling M[D] (Interval) = ", M_D)
        print("modeling M[N] (Interval) = ", M_N)
        lambdas_out_list.append(lambda_out)
        list_D.append(M_D)
        list_N.append(M_N)


    plt.figure(1)
    plt.xlabel('lambda')
    plt.axvline(lambd_critical, label='lambda critical = ' + str(round(lambd_critical,3)), color='r')
    plt.plot(lambdas, list_D, label ='M_D (Interval)')
    plt.legend()
    plt.show()
    plt.figure(2)
    plt.xlabel('lambda')
    plt.axvline(lambd_critical, label='lambda critical = ' + str(round(lambd_critical,3)), color='r')
    plt.plot(lambdas, list_N, label ='M_N (Interval)', color='darkorange')
    plt.legend()
    plt.show()
    plt.figure(3)
    plt.xlabel('lambda in')
    plt.ylabel('lambda out')
    plt.axhline(lambd_critical, label='lambda critical = ' + str(round(lambd_critical,3)), color='r')
    plt.plot(lambdas, lambdas_out_list, label='lambda out (Interval)', color = 'darkmagenta')
    plt.legend()
    plt.show()
    plt.close()
    
if __name__ == "__main__":
	main()
