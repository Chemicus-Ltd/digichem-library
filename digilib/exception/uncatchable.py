import signal

class Uncatchable_exception(BaseException):
    """
    Superclass for exceptions that are not normally caught.
    
    Note that no attempt is made to stop you catching these exceptions if you want, but they will not be caught by 'except Exception' clauses.
    """
    
class Submission_paused(Uncatchable_exception):
    """
    Exception raised during some submission routines.
    
    This exception signals the program that nothing further is to be done from this process; submission will continue in a second process (possible on a different machine on a different planet).
    
    You do not normally want to catch this exception; it is raised during the normal submission process.
    """
    
class Signal_caught(Uncatchable_exception):
    """
    Exception raised when a kill-type signal is sent and caught by the process.
    
    Similarly to KeyboardInterrupt; this exception is not normally caught (so that we exit NOW, this is important because some methods (SLURM, for example) time how long it takes us to shut down, and will send SIGKILL if we're too slow (which is uncatchable by design).
    However, it is valid (and expected) that vital cleanup handlers catch this exception, so long as they re-raise it once done.
    """
    
    def __init__(self, signalnum, stack_frame):
        """
        Constructor for Signal_caught exceptions.
        
        :param signalnum: The signal that was raised.
        :param stackframe: The current stack frame where the signal was caught (can be None).
        """
        self.signalnum = signalnum
        self.stack_frame = stack_frame
        
    @property
    def signal_name(self):
        """
        Get the string name (SIGINT etc) of the signal we caught.
        """
        return self.signal_to_name(self.signalnum)
        
    @classmethod
    def signal_to_name(self, signalnum):
        """
        Convert a singal number to a string (unrecognised signals will return "???"
        """
        try:
            return signal.Signals(signalnum).name
        except Exception:
            return "???"
        
    def __str__(self, *args, **kwargs):
        return "Received signal {} ({})".format(self.signalnum, self.signal_name)
        
        
    @classmethod
    def raise_from_signal(self, signalnum, stack_frame):
        """
        Method designed to be called by a signal handler; automatically raises a Signal_caught exception.
        """
        raise self(signalnum, stack_frame)