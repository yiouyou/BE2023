import sys
from loguru import logger
from bioelectronica._const import LOG_NORMDATA


logger.remove()

_format = "<green>[{extra[job]}]</green> <blue>{level}</blue> <green>{time:YYYY-MM-DD HH:mm:ss}</green> <yellow>{message}</yellow>"
logger.add(sys.stdout, level="INFO", format=_format, colorize=True)

logger.add(LOG_NORMDATA, level='DEBUG', format=_format, filter=lambda x: x["extra"]["job"] == "normdata")

logger_normdata = logger.bind(job="normdata")

