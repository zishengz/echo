# Initialize a flag to track if the message has been printed
_initialized = False

# Define a function to print the message once
def initialize_module():
    global _initialized
    if not _initialized:
        # logo art adapted from https://patorjk.com/software/taag/#p=testall&f=Graffiti&t=EChO
        print("""
 ███████        ██████ ██               ██████  
 ██            ██      ██              ██    ██ 
 █████         ██      ███████         ██    ██ 
 ██            ██      ██   ██         ██    ██ 
 ███████ lectro ██████ ██   ██ emical   ██████ ptimizer
        """)

        print('             Copyright © 2024 Zisheng Zhang')
        print('       Please cite: JACS, 2024, 146, 14, 9623-9630  ')
        print('     & JACS 2025, accepted. DOI:10.1021/jacs.5c00775')

        _initialized = True

# Call the function during the first import
initialize_module()

# The rest of your module initialization can go here if needed