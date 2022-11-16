#########################################################
#
# author: Yoni Leibner
# description: a Gui based on matplotlib and tkinter.
#               this tool meant to help peeling exponential
#               voltage responses
# date of modification: 16.11.2022
#
#########################################################

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import tkinter as Tkinter
import tkinter.filedialog as tkFileDialog
import os,pickle

SHIFT_FACTOR = 0.05 # up down shift in mV

class Peeler(object):
    def __init__(self):
        self.fig, self.ax = plt.subplots(1, 2)
        plt.subplots_adjust(bottom=0.2, wspace=0.35)
        self.buttoms_ax = plt.axes([0, 0.005, 0.1, 0.075])
        self.start_buttons()
        self.ready_to_fit = False
        plt.show()

    def start_buttons(self):
        b1 = plt.axes([0.0, 0.005, 0.2, 0.075])
        b2 = plt.axes([0.2, 0.005, 0.1, 0.075])
        b3 = plt.axes([0.3, 0.005, 0.1, 0.075])
        b4 = plt.axes([0.4, 0.005, 0.1, 0.075])
        b5 = plt.axes([0.5, 0.005, 0.1, 0.075])
        b6 = plt.axes([0.6, 0.005, 0.1, 0.075])
        b7 = plt.axes([0.7, 0.005, 0.1, 0.075])
        b8 = plt.axes([0.8, 0.005, 0.1, 0.075])
        b9 = plt.axes([0.9, 0.005, 0.1, 0.075])

        self.bopen = Button(b1, 'Open Saved')
        self.bopen_id = self.bopen.on_clicked(self.open_save)

        self.bload = Button(b2, 'Load')
        self.bload_id = self.bload.on_clicked(self.open)

        self.bcut = Button(b3, 'Cut')
        self.bcut_id = self.bcut.on_clicked(self.cut)

        self.bflip = Button(b4, 'Flip')
        self.bflip_id = self.bflip.on_clicked(self.flip)

        self.bup = Button(b5, 'Up')
        self.bup_id = self.bup.on_clicked(self.up)

        self.bdown = Button(b6, 'Down')
        self.bdown_id = self.bdown.on_clicked(self.down)

        self.bfit = Button(b7, 'peel')
        self.bfit_id = self.bfit.on_clicked(self.fit_btn)

        self.bsave = Button(b8, 'Save')
        self.bsave_id = self.bsave.on_clicked(self.save)

        self.breturn = Button(b9, 'return')
        self.breturn_id = self.breturn.on_clicked(self.return_)

    def set_buttons(self):
        self.bopen.label.set_text('Open Saved')
        self.bload.label.set_text('Load')
        self.bcut.label.set_text('Cut')
        self.bdown.on_clicked(self.cut)
        self.bflip.label.set_text('Flip')
        self.bdown.on_clicked(self.flip)
        self.bup.label.set_text('Up')
        self.bdown.on_clicked(self.up)
        self.bdown.label.set_text('Down')
        self.bdown.on_clicked(self.down)
        self.bfit.label.set_text('peel')
        self.bsave.label.set_text('Save')
        self.breturn.label.set_text('return')


    def set_buttons2(self):
        self.bopen.label.set_text('Open Saved')
        self.bload.label.set_text('Load')
        self.bcut.label.set_text('')
        self.bcut.on_clicked(self.pass_func)
        self.bflip.label.set_text('')
        self.bflip.on_clicked(self.pass_func)
        self.bup.label.set_text('')
        self.bup.on_clicked(self.pass_func)
        self.bdown.label.set_text('')
        self.bdown.on_clicked(self.pass_func)
        self.bfit.label.set_text('peel')
        self.bsave.label.set_text('Save')
        self.breturn.label.set_text('return')


    def pass_func(self, to_print):
        pass

    def open(self, event):
        self.V = []
        self.T = []
        self.tau = []
        self.C = []
        self.bounds = []
        self.action_names = ['load']

        Tkinter.Tk().withdraw()  # Close the root window
        in_path = tkFileDialog.askopenfilename()

        try:
            M=np.loadtxt(in_path)
            if M.shape[1] == 3:
                self.T.append(M[:,0])
                self.V.append(M[:,2] - M[:50,2].mean()) #removing V_rest

            elif M.shape[1] == 1:
                self.T.append(M[:, 0])
                self.V.append(M[:, 1] - M[:50, 1].mean())  # removing V_rest

            self.dt = self.T[0][1] - self.T[0][0]
            self.plot(0)
            self.ready_to_fit = True
            self.idx_fit_start = None
        except:
            print(in_path,'\nnot a valid file for peeling')


    def open_save(self, event):
        Tkinter.Tk().withdraw()  # Close the root window
        in_path = tkFileDialog.askopenfilename()

        try:
            data=pickle.load(open(in_path, 'rb'))
            self.V = data['V']
            self.T = data['T']
            self.tau = data['tau']
            self.C = data['C']
            self.bounds = data['bounds']
            self.action_names = data['action']
            self.idx_fit_start=data['idx_fit_start']
            self.ready_to_fit = data['ready_to_fit']
            current_num = len(self.V) - 1
            self.plot(current_num)
        except:
            print(in_path,'\nnot a valid file for peeling')


    def plot(self, num, color='k'):
        self.ax[0].clear()
        self.ax[1].clear()
        if self.fit_started():
            self.set_buttons2()
        else:
            self.set_buttons()
        self.ax[0].plot(self.T[num], self.V[num], color=color, zorder=1)
        self.ax[1].plot(self.T[num], np.log(self.V[num]), color=color, zorder=1)
        self.ax[0].set_xlabel('time (ms)')
        self.ax[0].set_ylabel('voltage (mV)')
        self.ax[1].set_xlabel('time (ms)')
        self.ax[1].set_ylabel('log-voltage (log(mV))')

        x_lim = self.ax[0].get_xlim()
        y_lim = self.ax[0].get_ylim()
        for i, tau in enumerate(self.tau):
            print(i, tau, y_lim[1]-(i+1)*(y_lim[1]-y_lim[0])/5.0)
            self.ax[0].text(x_lim[1]-(x_lim[1]-x_lim[0])/3.0, y_lim[1]-(i+1)*(y_lim[1]-y_lim[0])/20.0, 'tau_'+str(i)+' = '+str(round(tau, 3)), fontsize=5)

        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

    def fit_started(self):
        for action in self.action_names:
            if action .startswith('fit-'):
                return True
        return False

    def cut(self, event):
        if self.fit_started():
            self.fig.suptitle('cut is disable after fit')
            return
        self.fig.suptitle('choose 2 points for start-end')
        current_num = len(self.V)-1
        cut_fit = self.fig.ginput(2)
        start = cut_fit[0][0]
        end = cut_fit[1][0]
        print(start, end)
        if start < end:
            self.T.append(self.T[current_num][int(start/self.dt):int(end/self.dt)]-start)
            self.V.append(self.V[current_num][int(start/self.dt):int(end/self.dt)])
            self.action_names.append('cut')
            self.plot(current_num+1)
            self.fig.suptitle('')
        else:
            self.fig.suptitle('you must choose a valid points where delta X > 0!!!')

    def flip(self, event):
        if self.fit_started():
            self.fig.suptitle('cut is disable after fit')
            return
        current_num = len(self.action_names)-1
        self.T.append(self.T[current_num])
        self.V.append(-self.V[current_num])
        self.action_names.append('flip')
        self.plot(current_num+1)
        self.fig.suptitle('')

    def return_(self, event):
        current_num = len(self.V)-1
        if self.ready_to_fit:
            if self.action_names[-1].startswith('confirm-fit'):
                self.fig.suptitle('cant return from fit!!!')
                return
            self.T.pop(-1)
            self.V.pop(-1)
            self.plot(current_num-1)
        else:
            self.plot(current_num)
            self.bfit.label.set_text('peel')
            self.ready_to_fit=True
            self.C.pop(-1)
            self.tau.pop(-1)
            self.bounds.pop(-1)

        self.action_names.pop(-1)
        self.fig.suptitle('')

    def up(self, event):
        if self.fit_started():
            self.fig.suptitle('cut is disable after fit')
            return
        current_num = len(self.action_names) - 1
        self.T.append(self.T[current_num])
        self.V.append(self.V[current_num]+SHIFT_FACTOR)
        self.action_names.append('up')
        self.plot(current_num + 1)
        self.fig.suptitle('')

    def down(self, event):
        if self.fit_started():
            self.fig.suptitle('cut is disable after fit')
            return
        current_num = len(self.action_names) - 1
        self.T.append(self.T[current_num])
        self.V.append(self.V[current_num]-SHIFT_FACTOR)
        self.action_names.append('down')
        self.plot(current_num + 1)
        self.fig.suptitle('')

    @staticmethod
    def remove_exp(T_, V_, C_, tau_):
        to_remove = C_ * np.exp(-T_ / tau_)
        return T_.copy(), np.array(V_ - to_remove).copy()

    @staticmethod
    def find_exp(T, V):
        print(T.shape, V.shape)
        Pol = np.polyfit(T, np.log(V), 1)
        return np.exp(Pol[1]), -1.0 / Pol[0]  # this are C and tau

    def fit_btn(self, event):
        if self.ready_to_fit:
            self.fit(event)
        else:
            self.confirm_fit(event)

    def fit(self, event):
        self.fig.suptitle('choose 2 points for start-end')
        current_num = len(self.V)-1
        if self.idx_fit_start is None:
            self.idx_fit_start = current_num
        cut_fit = self.fig.ginput(2)
        start = cut_fit[0][0] #- self.T[current_num][0]
        end = cut_fit[1][0] #- self.T[current_num][0]
        if start < end:
            fit_idx = np.logical_and(self.T[current_num] >= start, self.T[current_num] <= end)
            C, tau = self.find_exp(self.T[current_num][fit_idx].copy(), self.V[current_num][fit_idx].copy())
            self.action_names.append('fit-'+str(len(self.C)))
            self.C.append(C)
            self.tau.append(tau)
            V_start = np.log(self.V[current_num][self.T[current_num]>=start][0])
            V_end = np.log(self.V[current_num][self.T[current_num]>=end][0])
            self.bounds.append([start, end])
            self.ax[1].scatter([start, end], [V_start, V_end], zorder=3, color='r')
            X = self.T[current_num].copy()
            self.ax[1].plot(X, np.log(self.C[-1]) + -1.0/self.tau[-1]*X, zorder=2)
            self.fig.suptitle('tau '+str(len(self.tau))+' = '+str(round(tau, 3)))

            self.temp_stuff = dict(start=start, end=end)
            self.ready_to_fit = False
            self.fig.canvas.draw()
            self.bfit.label.set_text('ok')
            # block other bottuns


    def confirm_fit(self, event):
        # unblock other bottuns
        start = self.temp_stuff['start']
        end = self.temp_stuff['end']
        current_num = len(self.V) - 1

        self.action_names.append('confirm-fit' + str(len(self.C)))
        T, V = self.remove_exp(self.T[current_num], self.V[current_num], self.C[-1], self.tau[-1])
        self.T.append(T)
        self.V.append(V)
        self.plot(current_num + 1)
        self.fig.suptitle('')
        self.bfit.label.set_text('peel')
        self.ready_to_fit = True

    def save(self, event):
        print('in save')
        if self.idx_fit_start is None:
            self.fig.suptitle('nothing to save, first you need to peel')
        print('in save: idx_fit_start=', self.idx_fit_start)
        base_idx = self.idx_fit_start
        self.plot(base_idx)
        print('in save', self.C, self.tau)
        X = self.T[base_idx].copy()
        y_lim = self.ax[1].get_ylim()
        estimated_V = None
        acc_c = 0
        for c, tau in zip(self.C, self.tau):
            acc_c+=c
            # self.ax[1].plot(X, c + -1.0 / tau * X, zorder=2)
            self.ax[1].plot(X, np.log(acc_c) + -1.0 / tau * X, zorder=2)
            if estimated_V is None:
                estimated_V  = c * np.exp(-X / tau)
            else:
                estimated_V += c * np.exp(-X/ tau)
        self.ax[0].plot(X, estimated_V, color='r', ls='--')
        self.ax[1].set_ylim(ymin=y_lim[0], ymax=y_lim[1])
        self.fig.suptitle('peeling results')
        self.fig.canvas.draw()
        f = tkFileDialog.asksaveasfile(mode='w', defaultextension=".png")
        if f:  # asksaveasfile return `None` if dialog closed with "cancel".
            name = f.name
            plt.savefig(name)
            f.close()
            data_save_name = os.path.join(os.path.dirname(name), '.'.join(os.path.basename(name).split('.')[:-1])+'_data.p')
            pickle.dump(dict(
                                T=self.T,
                                V=self.V,
                                action=self.action_names,
                                ready_to_fit=self.ready_to_fit,
                                bounds = self.bounds,
                                C=self.C,
                                tau=self.tau,
                                idx_fit_start=self.idx_fit_start,
                            ), open(data_save_name, 'wb'))
        else:
            self.fig.suptitle('didnt save, please try again')

a=Peeler()