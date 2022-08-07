import numpy as np
import scipy.stats as stats
from thresholdmodeling import thresh_modeling
import block_print


# if overlap:
# gen_dist is fitted over the whole range but used only below thresh. tail_model is both fitted and used above thresh
# else
# gen_dist and tail are used and fitted below and above thresh respectively. Not really needed - special case.
class SplitModel:
    def __init__(self, main_model, tail_model, thresh, tcoef, overlap=True, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.main_model = main_model
        self.tail_model = tail_model
        self.thresh = thresh
        if overlap:
            # c1 + c2 != 1 !!!!
            self.c1 = (1 - tcoef)/self.main_model.cdf(thresh)
            self.c2 = tcoef
        else:
            # c1 + c2 = 1
            self.c1 = 1 - tcoef
            self.c2 = tcoef
        self._ppfvec = np.vectorize(self._ppf_single, otypes='d')
        self._pdfvec = np.vectorize(self._pdf_single, otypes='d')

    def _pdf_single(self, x):
        if x < self.thresh:
            return self.c1*self.main_model.pdf(x)
        else:
            return self.c2*self.tail_model.pdf(x)

    def pdf(self, x):
        return self._pdfvec(x)

    def _ppf_single(self, q):
        if q < (1 - self.c2):
            return self.main_model.ppf(q/self.c1)
        else:
            return self.tail_model.ppf((q+self.c2-1)/self.c2)

    def ppf(self, q):
        return self._ppfvec(q)


def fit_dist(data, thresh, gen_dist=stats.invgauss, main_model=None, overlap=True):
    if main_model is None:
        main_data = data if overlap else np.delete(data, data > thresh)
        gparams = gen_dist.fit(main_data)
        main_model = gen_dist(*gparams)
    block_print.block_print()
    shape, scale, *_ = thresh_modeling.gpdfit(data, thresh, 'mle')
    block_print.enable_print()
    loc = thresh
    tail_model = stats.genpareto(shape, loc, scale)
    tail = np.delete(data, data < thresh)
    return SplitModel(main_model, tail_model, thresh, tail.size / data.size, overlap=overlap)


# general split model. Only pdf is provided.
# Models is a list of tuples (model, is_bounded). If is_bounded 100 % probability is assumed to be in respective range.
# Requires model.pdf. If not is_bounded also requires model.ppf.
# data (np.array) is used to get coefficiens
class GenSplitModel:
    def __init__(self, data, models, bounds):
        coefs = []
        bounds = [-float("inf")] + bounds + [float("inf")]
        for (model, is_bounded), left, right in zip(models, bounds, bounds[1:]):
            c = sum(np.logical_and(left < data, data <= right))/data.size
            if is_bounded:
                coefs.append(c)
            else:
                print(0)
                print(model.ppf.size)
                coefs.append(c/(model.ppf(right) - model.ppf(left)))
        self.models = [model for model, is_bounded in models]
        self.bounds = bounds
        self.coefs = coefs
        self._pdfvec = np.vectorize(self._pdf_single, otypes='d')

    def _pdf_single(self, x):
        for i in range(len(self.models)):
            if self.bounds[i] < x <= self.bounds[i+1]:
                return self.coefs[i]*self.models[i].pdf(x)

    def pdf(self, x):
        return self._pdfvec(x)
