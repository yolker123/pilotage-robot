class CalculatorLogic:
    def evaluate_expression(self, expression):
        try:
            return eval(expression)
        except Exception as e:
            return f"Erreur : {e}"
