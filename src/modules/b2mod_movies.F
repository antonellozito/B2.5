      module b2mod_movies
      use b2mod_types
      use b2mod_constants
#ifndef NO_CDF
      integer, save :: ncvdim, nfcdim, nsdim, natmdim, nmoldim, ntridim,
     & timedim, ntstep

      contains

      subroutine create_cdf_start(ncid, filename,
     & nCv, nFc, ns, natm, nmol, ntri)
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer ncid, nCv, nFc, ns, natm, nmol, ntri
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer ntstepid, timesid
* enter define mode
      iret = nf_create (filename, NF_WRITE, ncid)
      call check_cdf_status(iret)
* define dimensions
      if (index(trim(filename), 'b2movies') .gt. 0) then
        iret = nf_def_dim(ncid, 'ncv', nCv, ncvdim)
        call check_cdf_status(iret)
        iret = nf_def_dim(ncid, 'nfc', nFc, nfcdim)
        call check_cdf_status(iret)
        iret = nf_def_dim(ncid, 'ns', ns, nsdim)
        call check_cdf_status(iret)
        if(natm.gt.0) then
          iret = nf_def_dim(ncid, 'natm', natm, natmdim)
          call check_cdf_status(iret)
        endif
        if(nmol.gt.0) then
          iret = nf_def_dim(ncid, 'nmol', nmol, nmoldim)
          call check_cdf_status(iret)
        endif
      endif
      if (index(trim(filename), 'eirenemovies') .gt. 0) then
        iret = nf_def_dim(ncid, 'ntri', ntri, ntridim)
        call check_cdf_status(iret)
        if(natm.gt.0) then
          iret = nf_def_dim(ncid, 'natm', natm, natmdim)
          call check_cdf_status(iret)
        endif
        if(nmol.gt.0) then
          iret = nf_def_dim(ncid, 'nmol', nmol, nmoldim)
          call check_cdf_status(iret)
        endif
      endif
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, timedim)
      call check_cdf_status(iret)
* define variables
      dims(1) = 1
      iret = nf_def_var (ncid, 'ntstep', NF_INT, 0, dims, ntstepid)
      call check_cdf_status(iret)
      dims(1) = timedim
      iret = nf_def_var (ncid, 'times', NF_DOUBLE, 1, dims, timesid)
      call check_cdf_status(iret)

      return
      end subroutine create_cdf_start

      subroutine create_cdf_field2(ncid, name, stype, long_name, units)
      implicit none
#     include <netcdf.inc>
      character*(*) name, long_name, units
* file id
      integer ncid, stype
* variable shapes
      integer dims(2)
* dimensions and variables
      integer iret
      integer z_cvfct_id
* define variable
      select case (stype)
      case(1)
        dims(1) = ncvdim
      case(2)
        dims(1) = nfcdim
      case default
        write(*,*) 'Undefined "stype" in create_cdf_field2 ',stype
        call xerrab ('Undefined stype value')
      end select
      dims(2) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 2, dims, z_cvfct_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvfct_id, 'long_name',
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvfct_id, 'units',
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_field2

      subroutine create_cdf_field3(ncid, name, stype, itype, long_name, 
     &      units)
      implicit none
#     include <netcdf.inc>
      integer ncid, stype, itype
      character*(*) name, long_name, units
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer z_cvfcst_id
      external xerrab
* define variable
      select case (stype)
      case(1)
        dims(1) = ncvdim
      case(2)
        dims(1) = nfcdim
      case default
        write(*,*) 'Undefined "stype" in create_cdf_field3 ',stype
        call xerrab ('Undefined stype value')
      end select
      select case (itype)
      case(1)
        dims(2) = nsdim
      case(2)
        dims(2) = natmdim
      case(3)
        dims(2) = nmoldim
      case default
        write(*,*) 'Undefined "itype" in create_cdf_field3 ',itype
        call xerrab ('Undefined itype value')
      end select
      dims(3) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 3, dims, z_cvfcst_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvfcst_id, 'long_name',
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvfcst_id, 'units',
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_field3

      subroutine create_cdf_field_tri(ncid, name, itype, long_name,
     &      units)
      implicit none
#     include <netcdf.inc>
      character*(*) name, long_name, units
* file id
      integer ncid, itype
* variable shapes
      integer dims(3)
* dimensions and variables
      integer iret
      integer z_trit_id
* define variable
      dims(1) = ntridim
      select case (itype)
      case(1)
        dims(2) = natmdim
      case(2)
        dims(2) = nmoldim
      case default
        write(*,*) 'Undefined "itype" in create_cdf_field_tri ',itype
        call xerrab ('Undefined itype value')
      end select
      dims(3) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 3, dims, z_trit_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_trit_id, 'long_name',
     1 len_trim(long_name), trim(long_name))
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_trit_id, 'units',
     1 len_trim(units), trim(units))
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_field_tri

      subroutine create_cdf_b2movies_fields(ncid, natm, nmol, eirene,
     & switch)
*     ..defines the fields of b2movies.nc selected by the cdfmovie_*
*       switches (the EIRENE fields only if eirene); shared by the
*       code, its truncation of the file and the conversion by b2us,
*       so that the three always write the same layout
      use b2mod_switches, only : switches
      implicit none
      integer ncid, natm, nmol
      logical eirene
      type (switches), intent(in) :: switch
      if (switch%write_b25_state_variables.ge.1) then
        call create_cdf_field3(ncid, 'na', 1, 1,
     1   'Ion density', 'm^-3')
        call create_cdf_field2(ncid, 'ne', 1,
     1   'Electron density', 'm^-3')
        call create_cdf_field3(ncid, 'ua', 1, 1,
     1   'Ion parallel velocity', 'm s^-1')
        call create_cdf_field2(ncid, 'te', 1,
     1   'Electron temperature', 'eV')
        call create_cdf_field2(ncid, 'ti', 1,
     1   'Ion temperature', 'eV')
        call create_cdf_field2(ncid, 'po', 1,
     1   'Electric potential', 'V')
      endif
      if (switch%write_b25_fluxes.ge.1) then
        call create_cdf_field3(ncid, 'fnax', 2, 1,
     1   'Poloidal particle flux', 's^-1')
        call create_cdf_field3(ncid, 'fnay', 2, 1,
     1   'Radial particle flux', 's^-1')
        call create_cdf_field2(ncid, 'fhex', 2,
     1   'Poloidal electron heat flux', 'W')
        call create_cdf_field2(ncid, 'fhey', 2,
     1   'Radial electron heat flux', 'W')
        call create_cdf_field2(ncid, 'fhix', 2,
     1   'Poloidal ion heat flux', 'W')
        call create_cdf_field2(ncid, 'fhiy', 2,
     1   'Radial ion heat flux', 'W')
        call create_cdf_field2(ncid, 'fchx', 2,
     1   'Poloidal electric current', 'A')
        call create_cdf_field2(ncid, 'fchy', 2,
     1   'Radial electric current', 'A')
      endif
      if (switch%write_b25_sources.ge.1) then
        call create_cdf_field3(ncid, 'rsana', 1, 1,
     1   'Volumetric ionisation rate', 'm^-3 s^-1')
        call create_cdf_field3(ncid, 'rsahi', 1, 1,
     1   'Volumetric ionisation energy source/sink', 'W m^-3')
        call create_cdf_field2(ncid, 'rsahisum', 1,
     1   'Total volumetric ionisation energy source/sink',
     2   'W m^-3')
        call create_cdf_field3(ncid, 'rrana', 1, 1,
     1   'Volumetric recombination rate', 'm^-3 s^-1')
        call create_cdf_field3(ncid, 'rrahi', 1, 1,
     1   'Volumetric recombination energy source/sink',
     2   'W m^-3')
        call create_cdf_field2(ncid, 'rrahisum', 1,
     1   'Total volumetric recombination energy source/sink',
     2   'W m^-3')
        call create_cdf_field3(ncid, 'rcxna', 1, 1,
     1   'Volumetric charge exchange rate', 'm^-3 s^-1')
        call create_cdf_field3(ncid, 'rcxhi', 1, 1,
     1   'Volumetric charge exchange energy source/sink',
     2   'W m^-3')
        call create_cdf_field2(ncid, 'rcxhisum', 1,
     1   'Total volumetric charge exchange energy source/sink',
     2   'W m^-3')
        call create_cdf_field3(ncid, 'rqahe', 1, 1,
     1   'Volumetric electron cooling rate', 'W m^-3')
        call create_cdf_field2(ncid, 'rqahesum', 1,
     1   'Total volumetric electron cooling rate', 'W m^-3')
        call create_cdf_field3(ncid, 'rqrad', 1, 1,
     1   'Volumetric line radiation rate', 'W m^-3')
        call create_cdf_field2(ncid, 'rqradsum', 1,
     1   'Total volumetric line radiation rate', 'W m^-3')
      endif
      if (switch%write_b25_transport.ge.1) then
        call create_cdf_field3(ncid, 'dna', 1, 1,
     1   'Particle density-driven diffusivity', 'm^-2 s^-1')
        call create_cdf_field3(ncid, 'dpa', 1, 1,
     1   'Particle pressure-driven diffusivity', 'm^-2 s^-1')
        call create_cdf_field2(ncid, 'hci', 1,
     1   'Ion thermal anomalous diffusivity', 'm^-2 s^-1')
        call create_cdf_field2(ncid, 'hce', 1,
     1   'Electron thermal anomalous diffusivity', 'm^-2 s^-1')
        call create_cdf_field3(ncid, 'vlax', 1, 1,
     1   'Poloidal anomalous pinch velocity', 'm s^-1')
        call create_cdf_field3(ncid, 'vlay', 1, 1,
     1   'Radial anomalous pinch velocity', 'm s^-1')
        call create_cdf_field3(ncid, 'vsa', 1, 1,
     1   'Anomalous viscosity', 'm kg^-1 s^-1')
      endif
      if (eirene) then
        if (switch%write_eirene_state_variables.ge.1) then
          if(natm.gt.0) then
            call create_cdf_field3(ncid, 'dab2', 1, 2,
     1       'Atom density', 'm^-3')
            call create_cdf_field3(ncid, 'tab2', 1, 2,
     1       'Atom temperature', 'eV')
          endif
          if(nmol.gt.0) then
            call create_cdf_field3(ncid, 'dmb2', 1, 3,
     1       'Molecular density', 'm^-3')
            call create_cdf_field3(ncid, 'tmb2', 1, 3,
     1       'Molecule temperature', 'eV')
          endif
        endif
        if (switch%write_eirene_fluxes.ge.1) then
          if(natm.gt.0) then
            call create_cdf_field3(ncid, 'rfluxa', 1, 2,
     1       'Atomic radial flux density', 'm^-2 s^-1')
            call create_cdf_field3(ncid, 'pfluxa', 1, 2,
     1       'Atomic poloidal flux density', 'm^-2 s^-1')
            call create_cdf_field3(ncid, 'refluxa', 1, 2,
     1       'Atomic radial energy flux density', 'W m^-2')
            call create_cdf_field3(ncid, 'pefluxa', 1, 2,
     1       'Atomic poloidal energy flux density', 'W m^-2')
          endif
          if(nmol.gt.0) then
            call create_cdf_field3(ncid, 'rfluxm', 1, 3,
     1       'Molecular radial flux density', 'm^-2 s^-1')
            call create_cdf_field3(ncid, 'pfluxm', 1, 3,
     1       'Molecular poloidal flux density', 'm^-2 s^-1')
            call create_cdf_field3(ncid, 'refluxm', 1, 3,
     1       'Molecular radial energy flux density', 'W m^-2')
            call create_cdf_field3(ncid, 'pefluxm', 1, 3,
     1       'Molecular poloidal energy flux density', 'W m^-2')
          endif
        endif
        if (switch%write_eirene_sources.ge.1) then
          if(natm.gt.0) then
            call create_cdf_field3(ncid, 'eirene_papl_sna', 1, 1,
     1       'Particle source from atom-plasma interactions',
     2       'm^-3 s^-1')
            call create_cdf_field3(ncid, 'eirene_pppl_sna', 1, 1,
     1       'Particle source from recombinations',
     2       'm^-3 s^-1')
            call create_cdf_field2(ncid, 'eirene_eael_she', 1,
     1       'Electron energy source from atom-plasma interactions',
     2       'W m^-3')
            call create_cdf_field2(ncid, 'eirene_epel_she', 1,
     1       'Electron energy source from recombinations',
     2       'W m^-3')
            call create_cdf_field2(ncid, 'eirene_eapl_shi', 1,
     1       'Ion energy source from atom-plasma interactions',
     2       'W m^-3')
            call create_cdf_field2(ncid, 'eirene_eppl_shi', 1,
     1       'Ion energy source from recombinations',
     2       'W m^-3')
            call create_cdf_field3(ncid, 'eneutrad', 1, 2,
     1       'Volumetric radiation rate due to atoms', 'W m^-3')
            call create_cdf_field2(ncid, 'eneutradsum', 1,
     1       'Total volumetric radiation rate due to atoms',
     2       'W m^-3')
          endif
          if(nmol.gt.0) then
            call create_cdf_field3(ncid, 'eirene_pmpl_sna', 1, 1,
     1       'Particle source from molecule-plasma interactions',
     2       'm^-3 s^-1')
            call create_cdf_field2(ncid, 'eirene_emel_she', 1,
     1       'Electron energy source from molecule-plasma ' //
     2       'interactions', 'W m^-3')
            call create_cdf_field2(ncid, 'eirene_empl_shi', 1,
     1       'Ion energy source from molecule-plasma ' //
     2       'interactions', 'W m^-3')
            call create_cdf_field3(ncid, 'emolrad', 1, 3,
     1       'Volumetric radiation rate due to molecules',
     2       'W m^-3')
           call create_cdf_field2(ncid, 'emolradsum', 1,
     1       'Total volumetric radiation rate due to molecules',
     2       'W m^-3')
          endif
        endif
      endif
      return
      end subroutine create_cdf_b2movies_fields

      subroutine create_cdf_eirenemovies_fields(ncid, natm, nmol)
*     ..defines the fields of eirenemovies.nc
      implicit none
      integer ncid, natm, nmol
      if (natm.gt.0) then
        call create_cdf_field_tri(ncid, 'pdena', 1,
     1   'Atom particle density', 'm^-3')
        call create_cdf_field_tri(ncid, 'edena', 1,
     1   'Atom energy density', 'eV m^-3')
        call create_cdf_field_tri(ncid, 'vxdena', 1,
     1   'Atom x-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_field_tri(ncid, 'vydena', 1,
     1   'Atom y-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_field_tri(ncid, 'vzdena', 1,
     1   'Atom z-directed momentum density', 'kg s^-1 m^-2')
      endif
      if(nmol.gt.0) then
        call create_cdf_field_tri(ncid, 'pdenm', 2,
     1   'Molecular particle density', 'm^-3')
        call create_cdf_field_tri(ncid, 'edenm', 2,
     1   'Molecular energy density', 'eV m^-3')
        call create_cdf_field_tri(ncid, 'vxdenm', 2,
     1   'Molecular x-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_field_tri(ncid, 'vydenm', 2,
     1   'Molecular y-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_field_tri(ncid, 'vzdenm', 2,
     1   'Molecular z-directed momentum density', 'kg s^-1 m^-2')
      endif
      return
      end subroutine create_cdf_eirenemovies_fields

      subroutine create_cdf_end(ncid)
      implicit none
#     include <netcdf.inc>
* file id
      integer ncid
* dimensions and variables
      integer iret, varid, index(1)
* leave define mode
      iret = nf_enddef (ncid)
      call check_cdf_status(iret)
* set the initial timestep index (0)
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      ntstep=0
      write(*,*) 'create_cdf_end: ',ntstep
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf_end

      subroutine put_cdf_start(ncid, filename, time)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) filename
      real (kind=R8) :: time(1)
      integer iret
      integer varid
      integer start(1), count(1), index(1)
c
      iret=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      ntstep=ntstep+1
      write(*,*) 'put_cdf_start: ',trim(filename),ntstep
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'times',varid)
      call check_cdf_status(iret)
      start(1)=ntstep
      count(1)=1
      iret=nf_put_vara_double(ncid,varid,start,count,time)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_start

      subroutine put_cdf_field2(ncid, name, z_cvfct, nCvFc)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      integer nCvFc
      real (kind=R8) :: z_cvfct(nCvFc)
      integer iret
      integer varid
      integer start(2), count(2)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=ntstep
      count(1)=nCvFc
      count(2)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_cvfct)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_field2

      subroutine put_cdf_field3(ncid, name, z_cvfcst, nCvFc, ns)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      integer nCvFc, ns
      real (kind=R8) :: z_cvfcst(nCvFc, ns)
      integer iret
      integer varid
      integer start(3), count(3)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=ntstep
      count(1)=nCvFc
      count(2)=ns
      count(3)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_cvfcst)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_field3

      subroutine put_cdf_field_tri(ncid, name, z_trit, ntri, ns)
      implicit none
#     include <netcdf.inc>
      integer ncid
      character*(*) name
      integer ntri, ns
      real (kind=R8) :: z_trit(ntri, ns)
      integer iret
      integer varid
      integer start(3), count(3)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=1
      start(3)=ntstep
      count(1)=ntri
      count(2)=ns
      count(3)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_trit)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_field_tri

      subroutine put_cdf_end(ncid)
      implicit none
#     include <netcdf.inc>
      integer ncid
      integer iret
c
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf_end

      subroutine get_cdf_ntstep(filename, ntstep, ncid)
      implicit none
#     include <netcdf.inc>
      character*(*) filename
      integer ntstep
      integer iret
      integer ncid, varid
      integer index(1)
c
      iret=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      end subroutine get_cdf_ntstep

      subroutine get_cdf(filename, ncid, ntstep, time, name, z_cvt, nCv)
      implicit none
#     include <netcdf.inc>
      character*(*) filename, name
      integer ntstep, nCv
      real (kind=R8) :: time(1), z_cvt(1:nCv)
      integer iret
      integer ncid, varid
      integer start(2), count(2)
c
      iret=nf_open(filename,NF_NOWRITE,ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'times',varid)
      call check_cdf_status(iret)
      start(1)=ntstep
      count(1)=1
      iret=nf_get_vara_double(ncid,varid,start,count,time)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=ntstep
      count(1)=nCv
      count(2)=1
      iret=nf_get_vara_double(ncid,varid,start,count,z_cvt)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      end subroutine get_cdf

c OLD interface

      subroutine create_cdf(filename, ncid, nCv, name, long_name, units)
      implicit none
#     include <netcdf.inc>
      character*(*) filename, name, long_name, units
      integer nCv
      integer index(1), ntstep
* variable shapes
      integer dims(2)
* dimensions and variables
      integer iret
      integer ncid, ncvdim, timedim
      integer ntstepid, timesid, z_cvt_id
* enter define mode
      iret = nf_create (filename, NF_NOCLOBBER, ncid)
      call check_cdf_status(iret)
* define dimensions
      iret = nf_def_dim(ncid, 'ncv', nCv, ncvdim)
      call check_cdf_status(iret)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, timedim)
      call check_cdf_status(iret)
* define variables
      dims(1) = 1
      iret = nf_def_var (ncid, 'ntstep', NF_INT, 0, dims, ntstepid)
      call check_cdf_status(iret)
      dims(1) = timedim
      iret = nf_def_var (ncid, 'times', NF_DOUBLE, 1, dims, timesid)
      call check_cdf_status(iret)
      dims(1) = ncvdim
      dims(2) = timedim
      iret = nf_def_var (ncid, name, NF_DOUBLE, 2, dims, z_cvt_id)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvt_id, 'long_name',
     1 len_trim(long_name), long_name)
      call check_cdf_status(iret)
      iret = nf_put_att_text (ncid, z_cvt_id, 'units',
     1 len_trim(units), units)
      call check_cdf_status(iret)

* leave define mode
      iret = nf_enddef (ncid)
      call check_cdf_status(iret)
* set the initial timestep index (0)
      ntstep=0
      write(*,*) 'create_cdf: ',filename,ntstep
      iret=nf_put_var1_int(ncid,ntstepid,index,ntstep)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine create_cdf

      subroutine put_cdf(ncid, filename, time, name, z_cvt, nCv)
      implicit none
#     include <netcdf.inc>
      character*(*) filename, name
      integer nCv
      real (kind=R8) :: time(1), z_cvt(nCv)
      integer ntstep
      integer iret
      integer ncid, varid
      integer start(2), count(2), index(2)
c
      iret=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'ntstep',varid)
      call check_cdf_status(iret)
      index(1)=1
      iret=nf_get_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
      ntstep=ntstep+1
      write(*,*) 'put_cdf: ',filename,ntstep
      iret=nf_put_var1_int(ncid,varid,index,ntstep)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,'times',varid)
      call check_cdf_status(iret)
      start(1)=ntstep
      count(1)=1
      iret=nf_put_vara_double(ncid,varid,start,count,time)
      call check_cdf_status(iret)
c
      iret=nf_inq_varid(ncid,name,varid)
      call check_cdf_status(iret)
      start(1)=1
      start(2)=ntstep
      count(1)=nCv
      count(2)=1
      iret=nf_put_vara_double(ncid,varid,start,count,z_cvt)
      call check_cdf_status(iret)
      iret = nf_close (ncid)
      call check_cdf_status(iret)
      return
      end subroutine put_cdf
#endif

      end module b2mod_movies
