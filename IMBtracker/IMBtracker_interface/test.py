from PyQt5.QtWidgets import *
from PyQt5.QtGui import *
import random


class Splash(QWidget):
    def __init__(self):
        super().__init__()

        # CREATE THE TABLE
        self.table = QTableView(self)  # SELECTING THE VIEW
        self.table.setGeometry(0, 0, 575, 575)
        self.model = QStandardItemModel(self)  # SELECTING THE MODEL - FRAMEWORK THAT HANDLES QUERIES AND EDITS
        self.table.setModel(self.model)  # SETTING THE MODEL
        self.table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.populate()

        self.table.doubleClicked.connect(self.on_click)

    def on_click(self, signal):
        row = signal.row()  # RETRIEVES ROW OF CELL THAT WAS DOUBLE CLICKED
        column = signal.column()  # RETRIEVES COLUMN OF CELL THAT WAS DOUBLE CLICKED
        cell_dict = self.model.itemData(signal)  # RETURNS DICT VALUE OF SIGNAL
        cell_value = cell_dict.get(0)  # RETRIEVE VALUE FROM DICT

        index = signal.sibling(row, 0)
        index_dict = self.model.itemData(index)
        index_value = index_dict.get(0)
        print(
            'Row {}, Column {} clicked - value: {}\nColumn 1 contents: {}'.format(row, column, cell_value, index_value))

    def populate(self):
        # GENERATE A 4x10 GRID OF RANDOM NUMBERS.
        # VALUES WILL CONTAIN A LIST OF INT.
        # MODEL ONLY ACCEPTS STRINGS - MUST CONVERT.
        values = []
        for i in range(10):
            sub_values = []
            for i in range(4):
                value = random.randrange(1, 100)
                sub_values.append(value)
            values.append(sub_values)

        for value in values:
            row = []
            for item in value:
                cell = QStandardItem(str(item))
                row.append(cell)
            self.model.appendRow(row)

        self.show()


if __name__ == '__main__':
    import sys

    app = QApplication(sys.argv)
    ex = Splash()
    sys.exit(app.exec_())
