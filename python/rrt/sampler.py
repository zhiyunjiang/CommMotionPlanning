import numpy as np

from enum import Enum
class SamplerType(Enum):
    UNIRAND = 1
    DETERMINISTIC = 2

class Sampler():
    def __init__(self, sampler_type, dest_freq, sequence=None):
        self.stype = sampler_type
        self.destFreq = dest_freq

        self.seq = sequence
        if self.seq is None:
            self.seq = []

        self.region = None

        self.dests = None

        self.useHalfCount = False

    def setParams(self, pppi):
        self.region = pppi.region
        self.d = len(self.region)//2
        goal_region = pppi.goalRegion
        self.dests = goal_region.getPoints()

    def sample(self, count):
        d_count = count
        if self.useHalfCount:
            d_count = count//2

        new_sample = -1*np.ones(self.d)
        if (d_count%self.destFreq) == 0 and (self.stype == SamplerType.UNIRAND):
            new_sample = self._sampleDestination()
        elif self.stype == SamplerType.UNIRAND:
            new_sample = self._sampleUnirand()
        elif self.stype == SamplerType.DETERMINISTIC:
            new_sample = self.seq[count%len(self.seq)]

        self.seq.append(new_sample)
        return new_sample


    def _sampleUnirand(self):
        sample = np.zeros(self.d)
        for i in range(self.d):
            v_max = self.region[2*i +1]
            v_min = self.region[2*i]
            sample[i] = np.random.uniform(low = v_min, high = v_max)

        return sample

    def _sampleDestination(self):

        idx = np.random.choice([i for i in range(len(self.dests))])
        return self.dests[idx,:]
    

