# set definitions

add_compile_definitions(SQLITE_THREADSAFE=1)

add_library(sqlite OBJECT sqlite3.c)
