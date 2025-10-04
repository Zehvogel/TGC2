import json
from itertools import combinations


class AltSetupHandler():
    _sm_ref: dict[str, float]
    _variations: list[float]
    _alt_setup: dict[str, dict[str, float]]
    _mirror: bool
    _combinations: bool


    def __init__(self, js: str, mirror: bool = True, combinations: bool=True):
        conf = json.loads(js)
        self._sm_ref = conf["SM"]
        self._variations = conf["variations"]
        self._mirror = mirror
        self._combinations = combinations
        self._alt_setup = {}
        self._make_alt_configs()


    def make_name(self, parameter: str, var: float) -> str:
        # need to get rid of - to use as c++ identifier later
        alt_name = f"{parameter}_{'pos' if var > 0. else 'neg'}_{abs(var):.0e}".replace("-", "m")
        return alt_name


    def _add_1d_var(self, par: str, var: float) -> None:
        conf = self._sm_ref.copy()
        name = self.make_name(par, var)
        conf[par] += var
        self._alt_setup[name] = conf


    def _add_2d_var(self, par1: str, par2: str, var1: float, var2: float) -> None:
        conf = self._sm_ref.copy()
        name1 = self.make_name(par1, var1)
        name2 = self.make_name(par2, var2)
        name = f"{name1}_{name2}"
        conf[par1] += var1
        conf[par2] += var2
        self._alt_setup[name] = conf


    def _make_alt_configs(self):
        for var in self._variations:
            for par in self._sm_ref:
                # make 1 var up
                self._add_1d_var(par, var)
                if self._mirror:
                    # make 1 var down
                    self._add_1d_var(par, -var)
            # make 2 vars
            if self._combinations:
                for par1, par2 in combinations(self._sm_ref, r=2):
                    # ++
                    self._add_2d_var(par1, par2, var, var)
                    if self._mirror:
                        # +-
                        self._add_2d_var(par1, par2, var, -var)
                        # -+
                        self._add_2d_var(par1, par2, -var, var)
                        # --
                        self._add_2d_var(par1, par2, -var, -var)


    def get_alt_setup(self):
        # would probably be cleaner to return a deepcopy :/
        return self._alt_setup


    def get_variations(self):
        return self._variations


    # TODO: refactor
    # this one returns also the mirrored ones
    def get_variations_ext(self):
        if not self._mirror:
            return self.get_variations()
        res = {}
        for name, vars in self._variations.items():
            new_vars = [-1 * i for i in vars] + vars
            res[name] = new_vars
        return res


    @staticmethod
    def get_var_from_name_1d(name: str):
        parts = name.split("_")
        num = parts[-1]
        num = num.replace("m", "-")
        sign = -1. if parts[-2] == "neg" else 1.
        return sign * float(num)


    def get_pars(self) -> list[str]:
        return list(self._sm_ref.keys())


    # TODO also use this above maybe 
    # TODO super ugly refactor candidate
    def get_variations_nd(self):
        nd = len(self._sm_ref.keys())
        for var in self._variations:
            for i in range(nd):
                r = [0., 0., 0.]
                r[i] += var
                yield self.get_name_nd(r), r
                if self._mirror:
                    r = [0., 0., 0.]
                    r[i] -= var
                    yield self.get_name_nd(r), r
                for j in range(nd):
                    if j <= i:
                        continue
                    r = [0., 0., 0.]
                    r[i] += var
                    r[j] += var
                    yield self.get_name_nd(r), r
                    if self._mirror:
                        r = [0., 0., 0.]
                        r[i] += var
                        r[j] -= var
                        yield self.get_name_nd(r), r
                        r = [0., 0., 0.]
                        r[i] -= var
                        r[j] += var
                        yield self.get_name_nd(r), r
                        r = [0., 0., 0.]
                        r[i] -= var
                        r[j] -= var
                        yield self.get_name_nd(r), r


    def get_name_nd(self, vars: list[float]):
        names = []
        for i, v in enumerate(vars):
            if v != 0:
                par = list(self._sm_ref.keys())[i]
                name = self.make_name(par, v)
                names.append(name)
        return "_".join(names)


    def get_n_variations(self):
        return len(list(self.get_variations_nd()))
        # n_vars = len(self._variations)
        # if self._mirror:
        #     return (2 + 4) * n_vars
        # else:
        #     return (1 + 1) * n_vars
