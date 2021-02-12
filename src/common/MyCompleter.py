'''
---------------------------
 Licensing and Distribution
---------------------------

Program name: Pilgrim
Version     : 2021.1
License     : MIT/x11

Copyright (c) 2021, David Ferro Costas (david.ferro@usc.es) and
Antonio Fernandez Ramos (qf.ramos@usc.es)

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"),
to deal in the Software without restriction, including without limitation
the rights to use, copy, modify, merge, publish, distribute, sublicense,
and/or sell copies of the Software, and to permit persons to whom the Software
is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included
in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.
---------------------------

*----------------------------------*
| Module     :  common             |
| Sub-module :  MyCompleter        |
| Last Update:  2020/02/03 (Y/M/D) |
| Main Author:  David Ferro-Costas |
*----------------------------------*

This module contains the MyCompleter class
as well as a function to set it
'''

#=============================================#
import glob
#=============================================#

class MyCompleter(object):
      '''
      Custom completer
      '''
      def __init__(self, options=[], sort=False):
          self.options = options
          if sort: self.options.sort()
      def complete(self, text, state):
          # on first trigger, build possible matches
          if state == 0:
              # cache matches (entries that start with entered text)
              if text:
                  self.matches = [s for s in self.options if s and s.startswith(text)]
              # no text entered, all matches possible
              else:
                  self.matches = self.options[:]
              # No options given, so list files in folder
              if len(self.matches) == 0:
                  self.matches = (glob.glob(text+'*')+[None])
          # return match indexed by state
          try:
              return self.matches[state]
          except IndexError:
              return None

def set_completer(options=[]):
    import readline
    completer = MyCompleter(options)
    readline.set_completer(completer.complete)
    readline.parse_and_bind

