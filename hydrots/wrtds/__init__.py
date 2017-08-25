def distance(x0, x1, dimensions):
    #example of wrapped distance function
    delta = np.abs(x0 - x1)
    delta = np.where(delta > 0.5 * dimensions, delta - dimensions, delta)
    return np.sqrt((delta ** 2).sum(axis=-1))

def round_sig(x, sig=2):
    return round(x, sig-int(floor(log10(abs(x))))-1)

def datetime2year(dt): 
    year_part = dt - datetime(year=dt.year, month=1, day=1)
    year_length = datetime(year=dt.year+1, month=1, day=1) - datetime(year=dt.year, month=1, day=1)
    return dt.year + year_part/year_length

month_starts = ['2001-01','2001-02','2001-03','2001-04','2001-05','2001-06','2001-07','2001-08','2001-09','2001-10',
               '2001-11','2001-12']
decimal_months = np.asarray([datetime2year(pd.to_datetime(x)) for x in month_starts]) % 2001

str_months    = ['Jan 1','Feb 1','Mar 1','Apr 1','May 1','Jun 1','Jul 1','Aug 1', 'Sep 1', 'Oct 1', 'Nov 1', 'Dec 1']


#TODO: find out if ln(q+1) is an appropriate way to avoid taking ln(0)
# add a term for retransformation bias



class WRTDS():

    def __init__(self, c_in, q_in, t_in, t_step=1/365.25, q_space=None, h_q=2.0, h_t=7.0, h_s=0.5, annual=False):
        """

        Parameters:
        ===========
            c_in (array): Array containing concentration observations
            q_in (array): Array containing discharge observations
            t_in (array): Array containing time of observations in decimal years
            h_q discharge half-window width
            h_t trend half-window width
            h_s seasonal half-window width

        Consider:
        =========
            set h_q to 1
        """
        self.annual = annual
        self.c_in   = np.asarray(c_in)
        self.q_in   = np.asarray(q_in)
        self.t_in   = t_in
        self.t_y    = np.asarray(t_in.apply(datetime2year))

        # define discharge (q), trend (t), and seasonal (s) half-widths 
        self.h_q    = h_q
        self.h_t    = h_t
        self.h_s    = h_s
        self.n      = len(c_in) # number of observations

        # calculate q space
        if not q_space:
            self.qq = np.logspace(-1, np.log10( round_sig(self.q_in.max(),2) ), 10)


        self.t_step = t_step

        if not annual:
            start       = self.start = self.t_y.min()
            self.t_y   = self.t_y - start
            #stop        = self.t_in[-1]
            stop        = self.t_y.max()
            steps       = int((stop-0)/t_step)

            self.tt = np.linspace(0, stop, num=steps)

        if annual:
            self.start = 0
            self.tt = np.linspace(0, 1, num=365)
            #self.t_y = self.t_y % 1  #drop the year from t_y


#    def datetime2year(self.dt): 
#        year_part = dt - datetime(year=dt.year, month=1, day=1)
#        year_length = datetime(year=dt.year+1, month=1, day=1) - datetime(year=dt.year, month=1, day=1)
#        return dt.year + year_part/year_length


    def q_weights(self):
        """ Calculate weights for discharge observations.

        Calculate weights of discharge observations using a tricube weight fuction.
        Each row of the array corresponds to a given q, and the columns contain the weights
        for each observation.

        return:
            an array of all q weights. 
        """
        h_q = self.h_q  # half-widht for Q


        q, q_hat = np.meshgrid(self.q_in, self.qq)

        # calculate distance
        d_q = np.abs(np.log(q+1) - np.log(q_hat+1))

        # calculate weight
        w_q = np.where(d_q < h_q, (1-(d_q/h_q)**3)**3, 0)

        return w_q

    def t_weight(self,t):
        # return weight array for given t

        h_t = self.h_t
        h_s = self.h_s

        # calculate the trend distance for observations relative to a specific time
        #if not self.annual:
        d_t = np.abs(t - self.t_y) 

        if self.annual: # this tweak wraps the t_space so that 0=1
            #d_t = np.where(d_t > 0.5, np.abs(d_t - 1), d_t)
            d_t = 1
            #d_t = np.minimum(d_t, np.abs(self.t_y + t - 1)) 

        # calculate seasonal distance trend
        d_s = np.minimum( np.ceil(d_t) - d_t, d_t - np.floor(d_t)) 

        #return d_t , t

        # calculate trend and seasonal weights using a tricube weight function
        w_t = np.where(d_t < h_t, (1-(d_t/h_t)**3)**3, 0)
        w_s = np.where(d_s < h_s, (1-(d_s/h_s)**3)**3, 0)

        #import pdb; pdb.set_trace()
        return w_t * w_s

    def q_weight(self,q):
        # return weight array for given t

        h_q = self.h_q  # half-widht for Q

        # calculate distance
        d_q = np.abs(np.log(self.q_in+1) - np.log(q+1))

        # calculate weight
        w_q = np.where(d_q < h_q, (1-(d_q/h_q)**3)**3, 0)

        return w_q



        return w_t * w_s

    #def save(self, filename):
    #    with open(filename, 'wb+') as f:
    #        pickle.dump(self, f, protocol=pickle.HIGHEST_PROTOCOL)

    def save(self, filename):
        with open(filename, 'wb+') as f:
            dill.dump(self, f)

    def load(filename):
        f = open(filename, 'rb')
        return dill.load(f)

    month_starts = ['2001-01','2001-02','2001-03','2001-04','2001-05','2001-06','2001-07','2001-08','2001-09','2001-10',
               '2001-11','2001-11']

    decimal_months = np.asarray([datetime2year(pd.to_datetime(x)) for x in month_starts]) % 2001

    str_months    = ['Jan 1','Feb 1','Mar 1','Apr 1','May 1','Jun 1','Jul 1','Aug 1', 'Sep 1', 'Oct 1', 'Nov 1', 'Dec 1']


    def contourf(self,analyte): #XXX duplicate
        t_space = test.tt + self.start
        tt, qq = np.meshgrid(t_space, self.qq)

        if self.annual:
            fig, ax = plt.subplots(figsize=(6,3))
            ax.set_xticks(decimal_months)
            ax.set_xticklabels(str_months, rotation='45')

        if not self.annual:
            fig, ax = plt.subplots(figsize=(.7*self.t_y.max(),3))

        plt.contourf(tt,qq,self.c_out, aspect='auto', interpolation='bicubic', cmap=cm.OrRd)

        plt.semilogy()
        plt.xlabel('Time')
        plt.ylabel('Discharge (CFS)')
        cb1 = plt.colorbar()

        cb1.set_label(analyte)


    def contourf2(self,analyte): #XXX duplicate
        t_space = test.tt + self.start
        tt, qq = np.meshgrid(t_space, self.qq)

        #if self.annual:

        plt.contourf(tt,qq,self.c_out, aspect='auto', interpolation='bicubic', cmap=cm.OrRd)

        plt.semilogy()
        plt.xlabel('Time')
        plt.ylabel('Discharge (CFS)')
        cb1 = plt.colorbar()

        cb1.set_label(analyte)


    def contourf_log(self,analyte='SSC (mg/L)'):
        if self.annual:
            fig, ax = plt.subplots(figsize=(4,3))

        if not self.annual:
            fig, ax = plt.subplots(figsize=(.7*self.t_y.max(),3))

        t_space = test.tt + self.start
        tt, qq = np.meshgrid(t_space, self.qq)

        pcm = plt.contourf(tt,qq,self.c_out, aspect='auto', interpolation='bicubic', cmap=cm.OrRd,
                           norm=colors.LogNorm(vmin=1, vmax=10000), vmin=1, vmax=10000,
                           levels=[1,10,100,1000,10000])
                           #levels=[1,5,10,50,100,500,1000,5000,10000]) 
                            #norm=colors.LogNorm(vmin=self.c_out.min(), vmax=self.c_out.max()))
                          #l

        plt.semilogy()
        plt.xlabel('Time')
        plt.ylabel('Discharge (CFS)')
        #plt.title('Station ' + str(self.station)
        cb1 = plt.colorbar(pcm)

        cb1.set_label(analyte)

    def contourf_flux_log(self):
        if self.annual:
            fig, ax = plt.subplots(figsize=(6,3))

        if not self.annual:
            fig, ax = plt.subplots(figsize=(.7*self.t_y.max(),3))

        t_space = test.tt + self.start
        tt, qq = np.meshgrid(t_space, self.qq)
        flux = qq*self.c_out*0.028316846592*1000 #flux in g per m3
        pcm = plt.contourf(tt,qq,flux, aspect='auto', interpolation='bicubic', cmap=cm.OrRd,
                           norm=colors.LogNorm(vmin=1, vmax=10000))

        plt.semilogy()
        plt.xlabel('Time')
        plt.ylabel('Discharge (CFS)')
        #plt.title('Station ' + str(self.station)
        cb1 = plt.colorbar(pcm)

        cb1.set_label('SS Flux (mg/L)')

    def plot_concentration(self, analyte):
        fig, ax = plt.subplots(figsize=(.7*self.t_y.max(),3))
        
        t_space = test.tt + self.start
        
        plt.plot(self.hourly_concentration(), t_space)
        
        
        #plt.semilogy()
        plt.xlabel('Time')
        plt.ylabel('{} (mg/L)'.format(analyte))
        plt.title('Station {}'.format(str(self.station)))
        
    #def sample_concentration(self):
    #    self.concentration_m = np.exp(self.c_m) * self.smear
    #    return self.concentration_m
                  
    def calc(self, method='ols'):
        """ 
        Calculate WRTDS
        
        Parameters:
        -----------
        
        Returns:
        --------
        r : array_like
            Array containing 
        Aray of q vs t. 
        """
        self.c_out = np.empty([self.qq.size, self.tt.size])
        self.w_n   = np.empty([self.qq.size, self.tt.size])
        self.ln_c = np.empty([self.qq.size, self.tt.size])
        
        #XXX trying wiht normalize
        if method=='ridge':
            reg = linear_model.Ridge(alpha=4.7) #alpha was determined from SSC at station 201
        
        elif method=='ols':
            # run OLS regression
            reg = linear_model.LinearRegression(n_jobs=-1)
        
        else:
            print('invalid regression choice')
            
        n = self.c_in.size # number of discrete sample observations
        
        # define observations for regression
        # ln(c)= A*1 + D*t + C*ln(Q) + D*sin(2*pi*t) + E*cos(2*pi*t)
        x = np.ones((n,5))
        #x = np.zeros((n,5)) #XXX
        x[:,1] = np.log(self.q_in + 1)
        x[:,2] = self.t_y
        x[:,3] = np.sin(2 * np.pi * self.t_y)
        x[:,4] = np.cos(2 * np.pi * self.t_y)
        self.x = x
        
        # create a vector containing the dependent variable
        y = np.log(self.c_in) 
        # these vars are used in the model prediction
        #x_h =  = np.ones((n,5))
        ln_q = np.log(self.qq)
        cos_t = np.cos(2 * np.pi * self.tt)
        sin_t = np.sin(2 * np.pi * self.tt)
        t    = self.tt # just so the code reads pretty
        
        self.ln_q = ln_q
        self.cos_t = cos_t
        self.sin_t = sin_t
        
        w_q = self.q_weights() # should be a x by n
        
        for j in range(len(self.tt)): # change to time range
        # for a given time
        
            w_t = self.t_weight(t[j])
            
            for i in range(len(self.qq)):
                
                w = w_t * w_q[i,:]
                
                if w.sum() == 0: import pdb; pdb.set_trace()

                reg.fit(x, y, sample_weight=w)
                c = reg.coef_    
                
                ###XXX put exp term back in
                
                ##self.c_out[i,j] = np.exp(c[0] + c[1]*t[j] + c[2]*ln_q[i] + c[3]*sin_t[j] + c[4]*cos_t[j])
                
                ln_c_ij = c[0] + c[1]*ln_q[i] + c[2]*t[j] + c[3]*sin_t[j] + c[4]*cos_t[j]
                
                
                #if(np.exp(ln_c_ij)>1e7):
                #    return x,y,w
                    
                self.ln_c[i,j] = ln_c_ij
                self.w_n[i,j]  = w.sum()
        
        self.c_m = c_m = np.zeros(self.n)
        
        for i in range(self.n):
            #calculte weights
            w_t = self.t_weight(self.t_y[i])
            w_q = self.q_weight(self.q_in[i])
            w = w_t*w_q
            reg.fit(x, y, sample_weight=w)
            c_m[i] = reg.predict(x[i].reshape(1,-1))
            
        e_t = (y - c_m)
        
        # create a function for calculating residuals
        self.e = e_t # these are all the residuals
        self.res = (e_t**2).sum()
        self.r2  = 1 - self.res / (y - y.mean()).sum() 
        
        self.smear = (w*np.exp(e_t)).sum()/w.sum()
        self.c_out = np.exp(self.ln_c) * self.smear
        
    def load_continuous_discharge(self, q_c, t_c):
        '''Loads full record of daily (or hourly) of discharge observations
        
        '''
        #df = pd.read_csv()
    
        self.q_c = q_c
        self.t_c = t_c # continuous discharge record
        # convert contininous discharge record to decimal year
        self.t_cy = np.asarray(t_c.apply(datetime2year))
    
    def interp_continuous_discharge(self):
        '''Interpolate discharge to fixed grid, correcting for leap years
        
        '''
        step = 1/(365.25*24) # one decimal hour
        start = self.t_cy[0]
        stop  = self.t_cy[-1]
        #steps = (stop-start)/step
        #import pdb; pdb.set_trace()
        t_c_temp  = np.arange(start, stop, step) 
        self.q_c = np.interp(t_c_temp, self.t_cy, self.q_c) 
        self.t_cy = t_c_temp - start
        #import pdb; pdb.set_trace()

        
    def calc_continuous(self):
        '''Loads full record of daily (or hourly) of discharge observations
        
        '''
        n = self.c_in.size # number of discrete sample observations
        reg = linear_model.LinearRegression(n_jobs=-1)
        
        x = np.ones((n,5))
        #x = np.zeros((n,5)) #XXX
        x[:,1] = np.log(self.q_in + 1)
        x[:,2] = self.t_y
        x[:,3] = np.sin(2 * np.pi * self.t_y)
        x[:,4] = np.cos(2 * np.pi * self.t_y)
        y = np.log(self.c_in)
        
        x_c = np.ones((len(self.q_c),5))
        #x = np.zeros((n,5)) #XXX
        x_c[:,1] = np.log(self.q_c + 1)
        x_c[:,2] = self.t_cy
        x_c[:,3] = np.sin(2 * np.pi * self.t_cy)
        x_c[:,4] = np.cos(2 * np.pi * self.t_cy)
        
        self.c_cm = np.zeros_like(self.q_c)
        
        #df = pd.read_csv()
        for q,t, i in zip(self.q_c, self.t_cy, np.arange(len(self.t_cy))):
            #calculte weights
            w_t = self.t_weight(t)
            w_q = self.q_weight(q)
            #return w_t, w_q
            w = w_t*w_q
            reg.fit(x, y, sample_weight=w)
            #import pdb; pdb.set_trace()
            self.c_cm[i] = reg.predict(x_c[i].reshape(1,-1))

            
    
        

