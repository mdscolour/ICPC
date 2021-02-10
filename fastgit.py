# -*- coding: utf-8 -*-
import os
from datetime import datetime
stringtime = datetime.now().strftime("%m/%d/%Y, %H:%M:%S") # current date and time
os.system("git add --all")
os.system('git commit -m "fast commit %s"'%stringtime)
os.system("git push")
#os.system('mkdir code')

