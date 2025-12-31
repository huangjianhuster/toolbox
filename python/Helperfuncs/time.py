import time
from functools import wraps

def _format_duration(seconds: float) -> str:
    units = [
        ("day",   86400),
        ("hour",  3600),
        ("min",   60),
        ("sec",   1),
    ]
    for name, scale in units:
        if seconds >= scale:
            value = seconds / scale
            plural = 's' if value >= 2 else ''
            if value >= 10:
                return f"{value:.1f} {name}{plural}"
            else:
                return f"{value:.2f} {name}{plural}"
    return f"{seconds:.3f} secs"


def timeit(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.perf_counter()
        try:
            return func(*args, **kwargs)
        finally:
            elapsed = time.perf_counter() - start
            print(f"<func: {func.__name__}> took {_format_duration(elapsed)}")
    return wrapper


if __name__ == "__main__":

    @timeit
    def test():
        time.sleep(1)

    test()
