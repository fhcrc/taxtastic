import sqlite3

# TODO: inherit errors from other database engines

class OperationalError(sqlite3.OperationalError):
    pass

class IntegrityError(sqlite3.IntegrityError):
    pass
