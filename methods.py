import numpy as np

from pykalman import KalmanFilter


class AccOutlierFilter(object):
    """ Added class to simplify the removal of outliers from a
    time-series measurement. Uses mahalanobis distance and a third order
    Kalman filter. """

    def __init__(self, observations, sampling_period): 
        """ Initialize class.

        Parameters
        ----------
        observations : [n_timesteps]
            matrix of observations
        sampling period : float
            sampling period
        """

        self.sampling_period = sampling_period
        self.observations = observations

        self.n_dim_state = 3
        self.n_timesteps = len(observations)

        # System composed of position, velocity, and acceleration.
        transition_matrix = (
            np.array([[ 1. , sampling_period , sampling_period**2 ],
                      [ 0. , 1.              , sampling_period    ],
                      [ 0. , 0.              , 1.                 ]])
        )

        # Only position is observed
        observation_matrix = np.array([1., 0., 0.])

        p, res, _, _, _ = np.polyfit(
            sampling_period*np.arange(5),
            observations[:5], 2, full=True)

        initial_state_mean = p[::-1]
        initial_state_covariance = np.eye(self.n_dim_state)
        transition_covariance = np.eye(self.n_dim_state)
        observation_covariance = np.array([res/5])

        self.kf = KalmanFilter(
            transition_matrices=transition_matrix,
            observation_matrices=observation_matrix,
            initial_state_mean=initial_state_mean,
            initial_state_covariance=initial_state_covariance,
            transition_covariance=transition_covariance,
            observation_covariance=observation_covariance)


    def _get_max_dist(self, t):
        # return np.inf
        if t < 25: return np.inf

        else:
            return (4*self.likelihoods[:t].std()
                    + self.likelihoods[:t].mean())

    def run_filter(self):
        """ Run filter in an online fashion, popping outliers if they
        exceed threshold for mahalanobis distance """


        self.likelihoods = np.ma.zeros(self.n_timesteps)
        self.f_state_means = np.zeros((self.n_timesteps,
                                       self.n_dim_state))
        self.f_state_cov = np.zeros((self.n_timesteps,
                                               self.n_dim_state,
                                               self.n_dim_state))

        # kf._initialize_parameters()


        for t in range(self.n_timesteps - 1):
            if t == 0:
                self.f_state_means[t] = self.kf.initial_state_mean
                self.f_state_cov[t] = self.kf.initial_state_covariance

            self.f_state_means[t + 1], self.f_state_cov[t + 1] = (
                self.kf.filter_update(
                    self.f_state_means[t],
                    self.f_state_cov[t],
                    self.observations[t + 1],
                    max_mahalanobis_dist = self._get_max_dist(t)
                )
            )

            self.likelihoods[t+1] = self.kf.curr_mahalobis_dist

            # Mask the current value if an outlier was removed
            if self.likelihoods[t+1] > self._get_max_dist(t):
                self.likelihoods[t+1] = np.ma.masked

        self.outliers = np.array(self.likelihoods.mask)

        return self
    
    def recalibrate(self):
        """ Optimize the initial conditions, ignoring flagged
        outliers """

        obs = np.ma.array(self.observations)
        obs.mask = self.outliers
        self.kf.em(obs)

    def smooth(self):

        obs = np.ma.array(self.observations)
        obs.mask = self.outliers
        means, cov = self.kf.smooth(obs)

        return means





    

if __name__ == "__main__":
    sampling_period = 0.023

    y_test = np.array([
       100774.,   84304.,   73794.,   69194.,   55889.,   45772.,
        41533.,   36665.,   32828.,   29662.,   27836.,   26166.,
        25268.,   25629.,   25660.,   25815.,   26764.,   27496.,
        27826.,   29383.,   29270.,   30920.,   31085.,   31828.,
        32344.,   33406.,   35757.,   36727.,   37181.,   38006.,
        38893.,   39470.,   39542.,   39336.,   40429.,   40821.,
        40058.,   38594.,   39419.,   38202.,   37511.,   35520.,
        34148.,   33829.,   32797.,   30827.,   29920.,   28517.,
        28012.,   26877.,   26496.,   26155.,   25609.,   24784.,
        24412.,   24206.,   24629.,   24474.,   24742.,   25794.,
        26815.,   27341.,   28218.,   29136.,   29879.,   31240.,
        32539.,   34097.,   35283.,   35304.,   33076.,   38501.,
        39140.,   40863.,   41028.,   40986.,   42523.,   43616.,
        43379.,   43224.,   42276.,   43750.,   42337.,   42141.,
        42812.,   42203.,   41347.,   41038.,   40718.,   40089.,
        40130.,   40110.,   39821.,   39078.,   38944.,   38645.,
        38583.,   38243.,   38140.,   38903.,   38532.,   38697.,
        38975.,   39192.,   40347.,   40605.,   41358.,   42729.,
        42018.,   43028.,   43854.,   44266.,   45431.,   45576.,
        46514.,   47412.,   48031.,   48041.,   48484.,   49505.,
        49557.,   49248.,   50341.,   50485.,   50960.,   49495.,
        50093.,   50289.,   49650.,   50671.,   51795.,   52249.,
        51166.,   51795.,   51310.,   51403.,   51981.,   52569.,
        52465.,   52527.,   51929.,   52414.,   52465.,   52960.,
        53074.,   53486.,   54095.,   54703.,   54683.,   54868.,
        54920.,   54023.,   56075.,   56828.,   57839.,   57045.,
        57859.,   57870.,   58489.,   58839.,   59241.,   60613.,
        61077.,   60758.,   57488.,   61892.,   61871.,   61510.,
        63047.,   63872.,   64594.,   64099.,   63841.,   66296.,
        64852.,   65852.,   61376.,   57344.,   57767.,   60778.,
        65481.,   67874.,   61294.,   70916.,   72247.,   72567.,
        68792.,   71267.,   71133.,   71246.,   71525.,   70421.,
        71401.,   72618.,   73134.,   72866.,   72928.,   73144.,
        73546.,   73660.,   74217.,   74640.,   75032.,   75104.,
        75547.,   75269.,   75475.,   75599.,   76785.,   77393.,
        76125.,   77001.,   76919.,   77579.,   78002.,   77424.,
        79085.,   79724.,   80219.,   78157.,   80281.,   81055.,
        79838.,   81952.,   82643.,   82736.,   82901.,   82457.,
        83437.,   83499.,   83272.,   85851.,   86232.,   84902.,
        87109.,   85221.,   84407.,   85490.,   87243.,   87521.,
        88089.,   87985.,   88233.,   88274.,   88192.,   86768.,
        88656.,   89749.,   89553.,   88243.,   89450.,   88584.,
        89027.,   89574.,   89759.,   91348.,   91874.,   91203.,
        90801.,   90615.,   91482.,   90481.,   91245.,   93648.,
        93534.,   93318.,   92523.,   92255.,   94040.,   92843.,
        93431.,   94844.,   94493.,   94462.,   95556.,   94277.,
        94473.,   94174.,   96752.,   95721.,   93957.,   85221.,
        74248.,   63037.,   58509.,   59974.,   55982.,   72391.,
       100094.,  104333.,  105869.,  106447.,  103961.,  103559.,
       103652.,  103054.,  101496.,  100733.,   99021.,   99681.,
        98083.,   98588.])
    # sampling_period = 0.54912
    # y_test = (60 * np.exp(-0.02*x) * np.sin(2*np.pi*x/24.) + 1E-6*x**3 -
    #           1E-4*x**2 - 0.1*x + 60. + 2*np.random.randn(*x.shape))

    # y_test[200] = 1.

    y_test = y_test[10:]

    test = AccOutlierFilter(y_test, sampling_period)
    test.run_filter()
    test.recalibrate()
    test.run_filter()
    smooth = test.smooth()[:,0]
    

    import matplotlib.pyplot as plt
    x = np.arange(len(y_test))*sampling_period

    fig, axmatrix = plt.subplots(nrows=2, sharex=True)

    axmatrix[0].plot(x, test.observations)
    axmatrix[0].plot(x[test.outliers], test.observations[test.outliers],
                     'go')
    axmatrix[0].plot(x, test.f_state_means[:,0])
    axmatrix[0].plot(x, smooth)


    test.likelihoods.mask = np.ma.nomask
    axmatrix[1].plot(x, test.likelihoods)
    axmatrix[1].plot(x[test.outliers], test.likelihoods[test.outliers],
                     'go')

    plt.show()
    