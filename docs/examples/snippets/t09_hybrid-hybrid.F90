! the sensor on the data around the face, once per bias (the right bias reads the mirrored cells)
call seth_sensor_lcr(ug(i-2:i+3), rmu, SETH_SENSOR_LCR_PSIC, sp)
call seth_sensor_lcr(ug(i+3:i-2:-1), rmu, SETH_SENSOR_LCR_PSIC, sm)
if (sp == 1 .and. sm == 1) then        ! both smooth: the linear scheme, no weights
   call seth_linear5(fp(i-2:i+3), fm(i-2:i+3), rp, rm)
else
   call seth_teno5(fp(i-2:i+3), fm(i-2:i+3), SETH_TENO_EPS, SETH_TENO5_CT, rp, rm)
   nnonlinear = nnonlinear + 1
endif
