import numpy as np
import search_space

ss = search_space.Search_Space

space = ss(0.1,np.array([[0,1],[0,1]]),8,0)
space.addPoint((0.5,0.5))
space.checkCover()
space.addPoint((0,1),5)
