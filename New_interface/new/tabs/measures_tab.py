from PyQt5 import Qt, QtGui
from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout

from New_interface.EXAMPLE_organisation.oscilloscope_widget import OscilloscopeWidget
from New_interface.new.widgets.robot_widget import RobotWidget
from New_interface.new.widgets.visualisation_measures_widget import VisualisationMeasuresWidget

class MeasuresTab(QWidget):
        def __init__(self):
            super().__init__()
            self.init_ui()

        def init_ui(self):
            main_layout = QHBoxLayout()
            right_layout = QVBoxLayout()

            # Left part
            visualisation_measure_widget = VisualisationMeasuresWidget()
            main_layout.addWidget(visualisation_measure_widget, stretch=3)

            # Right part
            oscilloscope_widget = OscilloscopeWidget()
            right_layout.addWidget(oscilloscope_widget)
            robot_widget = RobotWidget(visualisation_measure_widget)
            right_layout.addWidget(robot_widget)
            right_layout.addStretch()


            main_layout.addLayout(right_layout, stretch=2)


            self.setLayout(main_layout)




