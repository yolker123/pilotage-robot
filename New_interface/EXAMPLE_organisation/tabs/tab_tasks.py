from PyQt5.QtWidgets import QWidget, QVBoxLayout, QPushButton, QListWidget, QLineEdit
from logic.task_logic import TaskLogic

class TasksTab(QWidget):
    def __init__(self):
        super().__init__()
        self.task_logic = TaskLogic()
        self.init_ui()

    def init_ui(self):
        layout = QVBoxLayout()
        self.task_list = QListWidget()
        self.task_input = QLineEdit()
        add_button = QPushButton("Ajouter Tâche")
        remove_button = QPushButton("Supprimer Sélection")

        layout.addWidget(self.task_input)
        layout.addWidget(add_button)
        layout.addWidget(remove_button)
        layout.addWidget(self.task_list)

        self.setLayout(layout)

        # Connexion des boutons
        add_button.clicked.connect(self.add_task)
        remove_button.clicked.connect(self.remove_task)

    def add_task(self):
        task = self.task_input.text()
        if task:
            self.task_logic.add_task(task)
            self.update_task_list()

    def remove_task(self):
        selected_items = self.task_list.selectedItems()
        for item in selected_items:
            self.task_logic.remove_task(item.text())
        self.update_task_list()

    def update_task_list(self):
        self.task_list.clear()
        self.task_list.addItems(self.task_logic.get_tasks())
