$(OBJDIR)/chkpoint.o             :src-usr/chkpoint.f;                          $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/IO_tools.o             :src-usr/IO_tools.f;                          $(F77) -c $(FL2) $< -o $@
$(OBJDIR)/uparam.o               :src-usr/uparam.f;                            $(F77) -c $(FL2) $< -o $@
