# -*- coding: utf-8 -*-

__author__ = "Benjamin Devillers (bdvllrs)"
__version__ = "1.0.5"

import os
import yaml

__all__ = ["config"]


def update_config(conf, new_conf):
    for item in new_conf.keys():
        if type(new_conf[item]) == dict and item in conf.keys():
            conf[item] = update_config(conf[item], new_conf[item])
        else:
            conf[item] = new_conf[item]
    return conf


def load_yaml(file):
    try:
        conf_dict = yaml.load(file, Loader=yaml.FullLoader)
    except AttributeError:
        conf_dict = yaml.load(file)
    return conf_dict


class Config:
    def __init__(self, path="config/", cfg_file="config_readcount.yaml", cfg=None, update=False):
        self.__is_none = False
        self.__data = cfg if cfg is not None else {}
        if (path is not None and cfg is None) or (cfg is not None and update is True):
            self.__path = os.path.abspath(os.path.join(os.curdir, path))
            with open(os.path.join(self.__path, cfg_file), "rb") as config_file:
                self.__data.update(load_yaml(config_file))

    def set_(self, key, value):
        self.__data[key] = value

    def set_subkey(self, key, subkey, value):
        self.__data[key][subkey] = value

    def values_(self):
        return self.__data

    def save_(self, file):
        file = os.path.abspath(os.path.join(os.curdir, file))
        with open(file, 'w') as f:
            yaml.dump(self.__data, f)

    def __getattr__(self, item):
        if type(self.__data[item]) == dict:
            return Config(cfg=self.__data[item])
        return self.__data[item]

    def __getitem__(self, item):
        return self.__data[item]


class Singleton:
    def __init__(self, cls):
        self.cls = cls
        self.instance = None

    def __call__(self, *args, **kwargs):
        if self.instance is None:
            self.instance = self.cls(*args, **kwargs)
            return self.instance

