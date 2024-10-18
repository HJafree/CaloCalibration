""" Mu2E Calorimeter Calibration
    Track event display version 2.3
    by Giacinto Boccia
    2024-09-03
    """

import ROOT as R
import numpy as np
from array import array
import crystalpos
#crystalpos.py is just a file with 2 np arrays with crystal x and y and the measure of the crystal side

class Hit:
    def __init__(self, number, x, y, time, v_val) -> None:
        self.n = number
        self.x = x
        self.y = y
        self.t = time
        self.v = v_val

class Crystal:
    side = crystalpos.crys_side
    
    def __init__(self, number : int, x : float, y : float) -> None:
        self.n = number
        self.x = x
        self.y = y
        #initializations
        self.hit_arr : list[Hit] = []
        self.hit_n : int = 0
        self.sum_v : float = 0.
        self.sum_t_diff :float = 0.

    def angles(self) -> tuple[np.array, np.array]:
        #Returs the angles of the crystal as [x_max, x_min], [y_max, y_min]
        x = np.zeros(2)
        y = np.zeros(2)
        x[0] = self.x + self.side / 2
        y[0] = self.y + self.side / 2
        x[1] = self.x - self.side / 2
        y[1] = self.y - self.side / 2
        return x, y
    
    def force_new_hit(self, new_hit : Hit) -> None:
        #hit_arr is only a temporary storage for sigle hits waiting for matching
        self.hit_arr.append(new_hit)
        if len(self.hit_arr) == 2:
            #If we reach 2 hits, the hit_arr is emptied to preserve memory
            self.hit_n += 2
            self.sum_v += self.hit_arr[1].v + self.hit_arr[0].v 
            t_diff = np.abs(self.hit_arr[1].t - self.hit_arr[0].t)
            if not np.isnan(t_diff):
                self.sum_t_diff += t_diff
            self.hit_arr = []
  
    def get_v(self) -> np.double:
        #An average is returned
        if self.hit_n > 0:
            return self.sum_v / self.hit_n
        else:
            return 0.
    
    def get_t_diff(self) ->np.double:
        #If more than two events, returns the mean of (|t1-t2|, |t3-t4|, |t5-t6|, ...)
        if self.hit_n > 0:
            return self.sum_t_diff / self.hit_n
        else:
            return np.nan

    def test_new_hit(self, new_hit : Hit) -> bool:
        #Check if hit is inside the crystal, get angles
        x, y = self.angles()
        x_min = np.min(x)
        x_max = np.max(x)
        y_min = np.min(y)
        y_max = np.max(y)
        if x_min <= new_hit.x <= x_max and y_min <= new_hit.y <= y_max:
            self.force_new_hit(new_hit)
            return True
        else:
            return False
        
    def empty(self) -> None:
        #Keeps crystal position and remove all hits and values
        self.hit_arr = []
        self.hit_n = 0
        self.sum_v = 0.
        self.sum_t_diff = 0.
        
        
class Disk:
    N_CRYSTALS = 674
       
    def __init__(self, id :int = 0) -> None:
        #cry_pos is a matrix where each row is the [x,y] position of a crystal
        self.fit_arr : list[Disk.Fit] =[]
        self.centroid : np.array[np.double, np.double] = np.zeros(shape = 2, dtype = np.double)
        
        #prepare the crystals
        self.cry_arr : list[Crystal] = []
        if id == 0:
            self.cry_pos = np.column_stack((crystalpos.crys_x, crystalpos.crys_y))[ : self.N_CRYSTALS-1]
        elif id == 1:
            self.cry_pos = np.column_stack((crystalpos.crys_x, crystalpos.crys_y))[self.N_CRYSTALS-1 : ]
        else:
            print("Invalid Disk id!")   
        for i in range(self.cry_pos.shape[0]):
            self.cry_arr.append(Crystal(i, self.cry_pos[i, 0], self.cry_pos[i, 1]))
        
    def load_event(self, slice) -> int:
        #slice is a root TTree slice, returns the number of loaded hits
        self.run_num : int = slice.nrun
        self.ev_num : int = slice.evnum
        hits_for_crystal : np.array[int] = np.zeros(shape = self.cry_pos.shape[0], dtype = int)
        n_hits : int = slice.nHits
        x_arr = np.frombuffer(slice.Xval)
        y_arr = np.frombuffer(slice.Yval)
        t_arr = np.frombuffer(slice.Tval) + np.frombuffer(slice.templTime)
        v_arr = np.frombuffer(slice.Vmax)
        cry_num_arr = np.frombuffer(slice.iCry, dtype= np.int32)
        
        #The centroid of the event is defined as the weighted average of the hit positions
        if n_hits > 0:
            self.centroid = np.array([np.average(x_arr, weights= v_arr), 
                                      np.average(y_arr, weights= v_arr)], 
                                     dtype= np.double)
        
        #Place hits in crystals
        for hit_num in range(n_hits):
            #Loop over hits
            curr_hit = Hit(hit_num, x_arr[hit_num], y_arr[hit_num], t_arr[hit_num], v_arr[hit_num])
            if self.cry_arr[cry_num_arr[hit_num] % self.N_CRYSTALS].test_new_hit(curr_hit):
                hits_for_crystal[cry_num_arr[hit_num] % self.N_CRYSTALS] += 1
            else:
                print("Error! Hit doesn't correspond to a crystal!")    
        #Check that each crystal has an even number of hits
        #If any of the crystals doesn't have an even number of hits, the mismathced hit is removed
        for crys_num in np.where(hits_for_crystal % 2 != 0)[0]:
            self.cry_arr[crys_num].hit_arr = []
                    
        return n_hits
    
    def event_fit(self, threshold : np.double = 4000. , type : str = 'linear') -> tuple[int, np.double]:
        #type indicates the kind of fit to perform, this function might be called more than once if more than one
        #fit is needed, all performed fits are stored in a fit_arr
        new_fit = self.Fit(self)
        self.fit_arr.append(new_fit)
        if type == 'linear':
            n_sel = new_fit.linear_fit(threshold)
            chi_squared = new_fit.chi_sqare
            return n_sel, chi_squared
        else:
            return None
        
    def draw_v(self, fits : bool = True, plot_name : str = False) -> None:
        #Use this method to draw the V values of an event, if "fits" it will display all the applied fits
        ev_name : str = plot_name or "Run " + str(self.run_num) + " Event " + str(self.ev_num)
        self.canvas = R.TCanvas(ev_name, ev_name, 1150, 1000)
        self.canvas.SetRightMargin(0.2)
        self.__draw_v(ev_name)
        self.__draw_circles()
        if fits:
            for fit in self.fit_arr:
                fit._Fit__fit_draw(ev_name)
        
    def __draw_circles(self) -> None:
        #Circes
        self.inner_c = R.TEllipse(0, 0, 374)
        self.inner_c.SetLineColor(2)
        self.inner_c.SetLineWidth(3)
        self.inner_c.SetFillStyle(0)
        self.inner_c.Draw()
        self.outer_c = R.TEllipse(0, 0, 660)
        self.outer_c.SetLineColor(2)
        self.outer_c.SetLineWidth(3)
        self.outer_c.SetFillStyle(0)
        self.outer_c.Draw()
    
    def __draw_v(self, ev_name : str) -> None:
        #Create an histogram for the crystals and write the V value of each crystal as Z
        self.crys_hist = R.TH2Poly()
        for crystal in self.cry_arr:
            boud_x, boud_y = crystal.angles()
            bin = self.crys_hist.AddBin(boud_x[0], boud_y[0], boud_x[1], boud_y[1])
            mean_v = crystal.get_v()
            if not np.isnan(mean_v):
                self.crys_hist.SetBinContent(bin, mean_v)
        #Draw the histogram as colored/empty boxes
        self.crys_hist.SetStats(0)
        self.crys_hist.SetTitle(ev_name)
        self.crys_hist.GetZaxis().SetTitle("Vmax")
        self.crys_hist.GetZaxis().SetTitleOffset(2.1)
        self.crys_hist.Draw('apl')
        self.crys_hist.Draw('zcol Cont0 same')
    
    def draw_tdif(self, plot_name : str = False) -> None:
        #Use this method to draw the mean time difference of the resposnses of each crystal 
        # (also over multiple events)
        name : str = plot_name or "Run " + str(self.run_num) + " Event " + str(self.ev_num) + " time differences"
        self.canvas = R.TCanvas(name, name, 1150, 1000)
        self.canvas.SetRightMargin(0.2)
        self.crys_hist = R.TH2Poly()
        for crystal in self.cry_arr:
            boud_x, boud_y = crystal.angles()
            bin = self.crys_hist.AddBin(boud_x[0], boud_y[0], boud_x[1], boud_y[1])
            t_diff = crystal.get_t_diff()
            if not np.isnan(t_diff):
                self.crys_hist.SetBinContent(bin, t_diff)
        #Draw the histogram as colored/empty boxes
        self.crys_hist.SetStats(0)
        self.crys_hist.SetTitle(name)
        self.crys_hist.GetZaxis().SetTitle("T Difference")
        self.crys_hist.GetZaxis().SetTitleOffset(2.1)
        self.crys_hist.Draw('apl')
        self.crys_hist.Draw('zcol Cont0 same')
        self.__draw_circles()
        
    def draw_hitcount(self) -> None:
        #use this method to draw a plot where each box represents the number of hits collected by each crystal 
        #(typically afther loading mutiple events)
        name : str = "Number of Hits"
        self.canvas = R.TCanvas(name, name, 1150, 1000)
        self.canvas.SetRightMargin(0.2)
        self.crys_hist = R.TH2Poly()
        for crystal in self.cry_arr:
            boud_x, boud_y = crystal.angles()
            bin = self.crys_hist.AddBin(boud_x[0], boud_y[0], boud_x[1], boud_y[1])
            self.crys_hist.SetBinContent(bin, crystal.hit_n)
        #Draw the histogram as colored/empty boxes
        self.crys_hist.SetStats(0)
        self.crys_hist.SetTitle(name)
        self.crys_hist.GetZaxis().SetTitle("# Hits")
        self.crys_hist.GetZaxis().SetTitleOffset(2.1)
        self.crys_hist.Draw('apl')
        self.crys_hist.Draw('zcol Cont0 same')
        self.__draw_circles()
        
    def empty(self) -> None:
        #Use this metod to empty the disk, deleting all the loaded hits, fits, etc. while keeping the crystals
        self.fit_arr : list[Disk.Fit] =[]
        self.ev_num : int = 0
        self.run_num : int = 0
        self.centroid : np.array[np.double, np.double] = np.zeros(shape = 2, dtype = np.double)
        self.hits_for_crystal : np.array[int] = np.zeros(shape = self.cry_pos.shape[0], dtype = int)
        for crystal in self.cry_arr:
            crystal.empty()
        
    class Fit:
        #Once you decide to fit an event, you instantiate this nested class, if new fit methods are developed, 
        #they should be placed in this class. To use them, call event_fit() on the Event, 
        #it returns an Event_fit object, if you vant to fit the same event vith different parametrs
        #call event_fit() more than once
        def __init__(self, disk : 'Disk') -> None:
            self.disk : 'Disk' = disk
            self.vertical : bool = False
            self.chi_sqare : np.double = np.nan

        def linear_fit(self, threshold : float) -> int:
            #The treshold applies on Vmax (mean value if crystal has 2 hits) 
            #and all hits over it get fitted linearly, returns the number of crystals considered
            n_sel = 0
            x_arr = array('f')
            y_arr = array('f')

            #Apply treshold
            for crystal in self.disk.cry_arr:
                if crystal.get_v() > threshold:
                    x_arr.append(crystal.x)
                    y_arr.append(crystal.y)
                    n_sel +=1

            if n_sel > 1:
                #Put points in a TGraph
                self.n_sel = n_sel
                err_arr = array ('f', np.full(n_sel, crystalpos.crys_side / 2))               
                self.graph = R.TGraphErrors(n_sel, x_arr, y_arr, err_arr, err_arr)
                
                #Check if the points are on a verticl line
                if not(self.__is_vertical(x_arr, y_arr)):
                    self.fit = R.TF1("Linear", "pol1")
                    self.graph.Fit("Linear")
                    self.chi_sqare = self.fit.GetChisquare()
                    #Set the color here, so that different fits (differnet methods) can have different colors
                    self.fit.SetLineColor(4)
                    self.fit.SetLineWidth(2)
                else:
                    center = self.disk.centroid
                    self.fit = R.TLine(center[0], 660, center[0], -660)
                    self.fit.SetLineColor(4)
                    self.fit.SetLineWidth(2)
            return n_sel

        def __is_vertical(self, x_arr : array, y_arr : array) -> bool:
            #This method is called when fitting (it is suposed to be shared between differnet fit tecniques),
            #if you want to access its result just read Event.Evnet_fit.vertical
            min_x = np.min(x_arr)
            max_x = np.max(x_arr)
            dx = max_x - min_x
            #If all xs are within a crystal side of the oters, the track is considerend vertical
            self.vertical = dx < crystalpos.crys_side
            return self.vertical

        def draw(self) -> None:
            #This method produces a TCanvas with the detector activation and the current fit 
            #(the one Event.Event_fit object on wich you call the method), if you are looking at a 
            #TCanvas with all the fits onverlayed over the calorimeter activation, 
            #or if you don't want to draw fits at all look at Event.event_draw()
            fit_name = "Run " + str(self.disk.run_num) + "Event " + str(self.disk.ev_num)
            self.canvas = R.TCanvas(fit_name, fit_name, 1150, 1000)
            self.canvas.SetRightMargin(0.2)
            self.disk.__draw_v(fit_name)
            self.__fit_draw(fit_name)
            self.disk.__draw_circles
                      
        def __fit_draw(self, name: str) -> None:
            #Please use draw() instead!
            if hasattr(self, 'graph'):
                #Style
                self.graph.SetMarkerColor(1)
                self.graph.SetMarkerSize(2)
                self.graph.GetXaxis().SetTitle('X (mm)')
                self.graph.GetYaxis().SetTitle('Y (mm)')
                self.graph.GetXaxis().SetLimits(-660, 660)
                self.graph.GetYaxis().SetRangeUser(-650, 650)
                self.graph.Draw('* same')
                if hasattr(self, 'fit'):
                    self.fit.Draw('same')
            else:
                print("You are drawing a fit that does not exist! Try linear_fit() [or other] before.")     

calo = Disk(0)

def single_event_v(tree : R.TTree,
                   run_n : int = 0,
                   ev_n : int = 0,
                   v_min :np.double = 300,
                   hits_min :int = 6,
                   chi_max : np.double = 20,
                   v_mode : str = 'i',
                   manual : bool = False) -> tuple[int, int]:
    
    def single_event_try(calo : Disk, tree : R.TTree, entry_num : int,
                         threshold : float = 4000,
                         min_hits : int = 6,
                         max_chi : float = 20,
                         v_mode :str = 'i') -> bool:
        #If the event matches the conditions is is shown and True is returned
        tree.GetEntry(entry_num)
        #Afther GerEntry, the tree objects gains all the attributes of a slice
        calo.load_event(slice = tree)
        hits, chi = calo.event_fit(threshold = threshold, type = 'linear')
        if v_mode == 'i':
            #Include vertical tracks if they meet max_chi
            if hits > min_hits and (chi < max_chi or calo.fit_arr[0].vertical):
                calo.draw_v()
                return True
            else:
                return False
        elif v_mode == 'e':
            #Exclude vertival tracks
            if hits > min_hits and chi < max_chi :
                calo.draw_v()
                return True
            else:
                return False
        else:
            #Only show vertical tracks
            if hits > min_hits and calo.fit_arr[0].vertical :
                calo.draw_v()
                return True
            else:
                return False
    
    global calo
    entry_n = tree.GetEntryNumberWithBestIndex(run_n, ev_n)
    if manual:
        #In this case whe just draw the selected event
        calo.empty()
        tree.GetEntry(entry_n)
        calo.load_event(slice = tree)
        calo.event_fit(threshold = v_min, type = 'linear')
        calo.draw_v()
    else:
        #If not in manual mode, we try untill an event meets our criteria
        while entry_n < tree.GetEntries() - 1:
            calo.empty()
            entry_n += 1
            if single_event_try(calo, tree, entry_n, v_min, hits_min, chi_max, v_mode):
                    break
                        
    tree.GetEntry(entry_n)
    return tree.nrun, tree.evnum

def signle_event_td(tree: R.TTree, run_num: int, event_num : int) -> None:
    #Display the reconstructed time difference for the channels of each crystal for a single event
    tree.GetEntryWithIndex(run_num, event_num)
    calo.load_event(tree)
    calo.draw_tdif()
    
def load_tree(tree : R.TTree) -> None:
    #Load the entire tree for average/cumulative drawings
    for entry_num, slice in enumerate(tree):
        calo.load_event(slice)
        if entry_num % 1000 == 0:
            print("Loading entry:", entry_num)

if __name__ == '__main__':
    #Parameters
    TRESHOLD = 4000
    HITS_MIN = 6
    CHI_MAX = 20

    #Open tree
    file_name = input("File to open:")
    file = R.TFile.Open(file_name)
    tree = file.sidet
    average_mode = input ("Display the average time diffeneces? y/[n]")
    hitcount_mode = input ("Display the total hit count? y/[n]")
    if average_mode == "y":
        load_tree(tree)
        calo.draw_tdif("Average Time Differences")
    elif hitcount_mode == "y":
        load_tree(tree)
        calo.draw_hitcount()
    else:
        entry_num : int = 0
        run_n : int = 0
        ev_n :int = 0
        while entry_num < tree.GetEntries():
            run_n, ev_n = single_event_v(tree, run_n, ev_n, v_min= TRESHOLD, hits_min= HITS_MIN, chi_max= CHI_MAX)
            input("Press any key for T Differences")
            signle_event_td(tree, run_n, ev_n)
            input("Press any key for next")     
    input("Press any key to exit")