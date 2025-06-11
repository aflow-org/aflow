#!/usr/bin/env python

import json
import subprocess
import os

VERBOSE=False   #CO20200520

class CHull:

    def __init__(self, aflow_executable = 'aflow'):
        self.aflow_executable = aflow_executable

    def aflow_command(self, cmd):
        try:
            if VERBOSE: print('aflow_command(): cmd = ' + self.aflow_executable + cmd)
            output = subprocess.check_output(   #MB20200409
                    self.aflow_executable + cmd,
                    shell=True
                    )
            if VERBOSE: print('aflow_command(): output = ' + output)
            return output
        except subprocess.CalledProcessError:
            raise AssertionError('aflow executable not found: ' + self.aflow_executable)

    def REST_API_down(self, cmd):   #CO20200520
        cmd = cmd.replace(" --screen_only","")
        if VERBOSE: print('REST_API_down(): cmd = ' + self.aflow_executable + cmd)
        try:
            output = subprocess.check_output(   #MB20200409
                    self.aflow_executable + cmd,
                    shell=True
                    )
            if VERBOSE: print('REST_API_down(): output = ' + output)
            if "REST-API appears to be down" in output.decode('utf-8'): #MB20200409
                return True
            return False
        except subprocess.CalledProcessError:
            raise AssertionError('aflow executable not found: ' + self.aflow_executable)

    def get_hull(self, alloy, options = None):
        command = ' --chull --force'
        if options:
            command += ' ' + options
        command += ' --print=json --screen_only --alloy=' + alloy

        output = self.aflow_command(command)
        if VERBOSE: print('get_hull(): output = ' + output)

        res_json = json.loads(output)

        if not res_json:
            if self.REST_API_down(command):
                raise AssertionError('AFLOW REST-API appears to be down')

        return res_json

    def get_distance_to_hull(self, alloy, off_hull_point, options = None):
        command = ' --chull --force --distance_to_hull=' + off_hull_point
        if options:
            command += ' ' + options
        command += ' --print=json --screen_only --alloy=' + alloy

        output = self.aflow_command(command)
        if VERBOSE: print('get_distance_to_hull(): output = ' + output)

        res_json = json.loads(output)

        if not res_json:
            if self.REST_API_down(command):
                raise AssertionError('AFLOW REST-API appears to be down')

        return res_json

    def get_stability_criterion(self, alloy, hull_point, options = None):
        command = ' --chull --force --stability_criterion=' + hull_point
        if options:
            command += ' ' + options
        command += ' --print=json --screen_only --alloy=' + alloy

        output = self.aflow_command(command)
        if VERBOSE: print('get_stability_criterion(): output = ' + output)

        res_json = json.loads(output)

        if not res_json:
            if self.REST_API_down(command):
                raise AssertionError('AFLOW REST-API appears to be down')

        return res_json

    def get_hull_energy(self, alloy, composition, options = None):
        command = ' --chull --force --hull_energy=' + ','.join([ str(comp) for comp in composition ])
        if options:
            command += ' ' + options
        command += ' --print=json --screen_only --alloy=' + alloy

        output = self.aflow_command(command)
        if VERBOSE: print('get_hull_energy(): output = ' + output)

        res_json = json.loads(output)

        if not res_json:
            if self.REST_API_down(command):
                raise AssertionError('AFLOW REST-API appears to be down')

        return res_json
