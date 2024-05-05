import traceback

def error(message):
    """ Print error message and exit program. """

    print('E R R O R: '+str(message))
        # spaces ' ' to avoid popup window in ChimeraX
    traceback.print_exc()
    raise Exception()

    return


def e(message):
    """
    Print information about a error.
    This function is intended to use before execution of error().
    """
    print(message)

    return


def warning(message):
    """ Print warning message. """

    print('WARNING: '+str(message))
    raise Error('WARNING: '+str(message))

    return


def notice(message):
    """ Print notice message. """

    print('NOTICE: '+str(message))
    raise Exception('NOTICE: '+str(message))

    return


def typecheck(var, vartype):
    """ Typecheck of given variable. """

    if not(type(var) is vartype):
        print('1)')
        print(type(var))
        print('2)')
        print(type(vartype))
        error('typecheck failed: '+str(var)+'1) is not 2)')

    return


def p(message):
    """ Print message. """

    print(message)

    return


def d(message):
    """ Print messages if debug mode enabled. """
    debug_mode = True # False

    if debug_mode:
        print(message)

    return
