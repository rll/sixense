import roslib
roslib.load_manifest("sixense")
import rospy
from sixense.msg import Calib

class HydraPub:
    def __init__(self):
        self.pub = rospy.Publisher("hydra_calib", Calib)
    def publish(self, msg):
        self.pub.publish(msg)
