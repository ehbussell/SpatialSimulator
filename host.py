class Host(object):

    def __init__(self, x, y, state=None, reg=0):
        self.x = x
        self.y = y
        self.state = state
        self.reg = reg

        self.jump_times = {str(state): 0.0}

    def __repr__(self):
        repr_str = "Host(" + str(self.x) + ", " + str(self.y) + ", '"
        repr_str += str(self.state) + "', " + str(self.reg) + ")"

        return repr_str
