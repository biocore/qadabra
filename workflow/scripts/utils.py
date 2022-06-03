import logging


def get_logger(logfile, rulename):
    qadabra_logger = logging.getLogger("qadabra")
    qadabra_logger.setLevel(logging.INFO)
    fh = logging.FileHandler(logfile, mode="w")
    formatter = logging.Formatter(
        f"[%(asctime)s - {rulename}] :: %(message)s"
    )
    fh.setFormatter(formatter)
    qadabra_logger.addHandler(fh)

    logging.captureWarnings(True)
    logging.getLogger("py.warnings").addHandler(fh)

    return qadabra_logger
