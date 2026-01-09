import time


class TimerError(Exception):
    """A custom exception used to report errors in use of Timer class"""


class pytimer:
    """timer class"""

    def __init__(self):
        self._start_time = None

    def start(self):
        """Start a new timer"""
        if self._start_time is not None:
            raise TimerError(f"Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()

    def stop(self):
        """Stop the timer, and report the elapsed time"""
        if self._start_time is None:
            raise TimerError(f"Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._start_time = None
        print(f"Elapsed time: {elapsed_time:0.4f} seconds")

    def now(self, iFlag_with_seconds=1):
        """Get the current time as a string"""
        t_struct = time.localtime()
        if iFlag_with_seconds == 1:
            sTime = time.strftime("%Y%m%d_%H%M%S", t_struct)
        else:
            sTime = time.strftime("%Y%m%d_%H%M", t_struct)
        return sTime
