"""
Name: Recoder

Coder: HaoLing ZHANG (BGI-Research)[V1]

Current Version: 1

Functions:
(1) Record the whole state of fusing two Pattern DNA.
"""

import copy

import numpy


# noinspection PyUnusedLocal
class Recoder(object):

    def __init__(self, major_fork_indices, minor_fork_indices, start_position):
        """
        Initialize recoder for save process variables in the process of calculating maximum overlap.

        :param major_fork_indices: index list with fork of major strand.
        :param minor_fork_indices: index list with fork of minor strand.
        :param start_position: the start position of major strand.
        """
        self._major_fork_indices = copy.deepcopy(major_fork_indices)
        self._minor_fork_indices = copy.deepcopy(minor_fork_indices)
        self._major_state = 0
        self._minor_state = 0

        self._start_position = start_position

        if start_position in major_fork_indices:
            self._major_state = 3

        if 0 in minor_fork_indices:
            self._minor_state = 3

        self._fork_chooser = numpy.array([[3 for index in range(int(len(major_fork_indices) / 3))],
                                          [3 for index in range(int(len(minor_fork_indices) / 3))]])

    def update(self, current_position, passes):
        """
        Update the current major and minor state by the current position and current passes.

        :param current_position: current reading position in minor strand/
        :param passes: current passes for two strand.
        """
        if self._major_state > 0:
            self._major_state = 2 * int(passes[1][0] or passes[1][1]) + int(passes[0][0] or passes[0][1])
        if self._minor_state > 0:
            self._minor_state = 2 * int(passes[0][1] or passes[1][1]) + int(passes[0][0] or passes[1][0])

        fork_situations = [
            numpy.where(self._major_fork_indices == self._start_position + current_position),
            numpy.where(self._major_fork_indices == self._start_position + current_position + 1),
            numpy.where(self._minor_fork_indices == current_position),
            numpy.where(self._minor_fork_indices == current_position + 1)
        ]

        if len(fork_situations[0]) > 0:
            fork_position = fork_situations[0][0]
            if fork_position % 3 == 2:
                self._fork_chooser[0][int(fork_position / 3)] = self._major_state
                self._major_state = 0
        if len(fork_situations[1]) > 0:
            fork_position = fork_situations[1][0]
            if fork_position % 3 == 0:
                self._major_state = 3

        if len(fork_situations[2]) > 0:
            fork_position = fork_situations[2][0]
            if fork_position % 3 == 2:
                self._fork_chooser[1][int(fork_position / 3)] = self._minor_state
                self._minor_state = 0
        if len(fork_situations[3]) > 0:
            fork_position = fork_situations[3][0]
            if fork_position % 3 == 0:
                self._minor_state = 3

    def get_current_state(self):
        """
        Get current major state and minor state.

        :return: current major state and minor state.
        """
        return self._major_state, self._minor_state

    def get_start_position(self):
        """
        Get the start position of this recoder.

        :return: start position.
        """
        return self._start_position

    def get_fuse_chooser(self):
        """
        Get the fuse chooser.

        :return: fuse chooser of major fork and minor fork.
        """
        return self._fork_chooser

    def __str__(self):
        """

        :return:
        """
        return "start position = " + str(self._start_position) + "\n" + \
               "major chooser = " + str(self._fork_chooser[0]) + "\n" + \
               "minor chooser = " + str(self._fork_chooser[1])
