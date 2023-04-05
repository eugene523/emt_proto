class ValidationResult:
    def __init__(self):
        self.issues = []

    def add_issue(self, issue_msg: str):
        self.issues.append(issue_msg)

    def is_ok(self):
        return len(self.issues) == 0
    
    def __str__(self):
        if self.is_ok():
            return "no issues"
        
        info = ""
        for i in self.issues:
            info += str(i) + '\n'
        return info