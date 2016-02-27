Module general
! --------------------------------------------------------------!
! Ben Palmer, University of Birmingham
! --------------------------------------------------------------!
  Use kinds
  Use strings
! Force declaration of all variables
  Implicit None
! Public variables
! Make private
  Private
! Public
! ---- Variables
! ---- Functions
  Public :: GetClockTime
  Public :: dpToString
  Public :: intToString
! ---- Subroutines
  Public :: makeDir
  Public :: rmFile
  Public :: rmDir
  Public :: correctFilePath
  Public :: readFile
  Public :: extractArrayColumnDP
  Public :: extractArrayColumnInt
  Public :: swapArrayRows1D
  Public :: swapArrayRows2D
  Public :: strToIntArr
  Public :: strToDPArr
  Public :: strToStrArr
  Public :: timeAcc
! Functions
  Public :: FileExists
  Public :: CountRowFields

! ---- Subroutines
!---------------------------------------------------------------------------------------------------------------------------------------
  Contains
!---------------------------------------------------------------------------------------------------------------------------------------


! ---------------------------------------------------------
! MODULE FUNCTIONS
! ---------------------------------------------------------

! TIMES

  Function GetClockTime () RESULT (outputTime)
! Force declaration of all variables
    Implicit None
    Real(kind=DoubleReal) :: outputTime
    Call cpu_time(outputTime)
  End Function GetClockTime


! TYPE CONVERSION

  Function dpToString(inputDP) RESULT (outputString)
! Force declaration of all variables
    Implicit None
! declare private variables
    Real(kind=DoubleReal) :: inputDP
    Character(len=32) :: outputString
! Read dp to string
    inputDP = 1.0D0 * inputDP
    Write(outputString,"(ES16.8E3)") inputDP
  End Function dpToString

  Function intToString(inputInt) RESULT (outputString)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: inputInt
    Character(len=32) :: outputString
! Read int to string
    Write(outputString,"(I16)") inputInt
  End Function intToString

!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------



! ---------------------------------------------------------
! MODULE SUBROUTINES
! ---------------------------------------------------------

! FILES AND DIRECTORIES

  Subroutine makeDir(path)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("mkdir -p "//trim(path))
  End Subroutine makeDir

  Subroutine rmFile(path)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("rm -f "//trim(path))
  End Subroutine rmFile

  Subroutine rmDir(path)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Character(*) :: path
    Call system("rm -fR "//trim(path))
  End Subroutine rmDir


  Subroutine correctFilePath (filePath)
! Correct file path
    Implicit None  !Force declaration of all variables
! Declare variables
    Character(*), Intent(INOUT) :: filePath
    Character(Len(filePath)) :: tempFilePath
    Integer(kind=StandardInteger) :: i, n
! Blank string
    tempFilePath = BlankString(tempFilePath)
    n = 0
    Do i=1,Len(filePath)
      If((iachar(filePath(i:i)).ge.32.and.iachar(filePath(i:i)).le.126))Then
        If(filePath(i:i).ne."?".and.filePath(i:i).ne."*".and.&
          filePath(i:i).ne."%".and.filePath(i:i).ne."+")Then
          n = n + 1
          tempFilePath(n:n) = filePath(i:i)
        End If
      End If
    End Do
    filePath = tempFilePath
  End Subroutine correctFilePath


! READING FILES

  Subroutine readFile(inputFilePath, fileArray, n)
! Subroutine to read file into an array
! Removes comments !.....
! Removes blank lines
! Removes leading spaces
    Implicit None ! Force declaration of all variables
! Private variables
    Character(*) :: inputFilePath
    Integer(kind=StandardInteger) :: i, j, n, ios
    Character(*), Dimension(:) :: fileArray
    Character(len=255) :: fileRow, fileRowTemp
    Logical :: commentsFlag
! Exit if file does not exist
    If(.not.FileExists(inputFilePath))Then
      print *,"File Does Not Exist"
      stop
    End If
! open file
    Open(UNIT=9999,FILE=trim(inputFilePath),status='old',action='read')
    n = 0
    Do i=1,size(fileArray,1)
      Read(9999,"(A255)",IOSTAT=ios) fileRow
      If(ios /= 0)Then
        EXIT
      End If
! remove comments
      commentsFlag = .false.
      Do j=1,255
        If(fileRow(j:j).eq."!")Then
          commentsFlag = .true.
        End If
        If(commentsFlag)Then
          fileRow(j:j) = " "
        End If
      End Do
! remove blank lines
      fileRowTemp = trim(adjustl(fileRow))
      If(fileRowTemp(1:1).ne." ")Then
        n = n + 1
        fileArray(n) = fileRowTemp
      End If
    End Do
! Close file
    close(9999)
  End Subroutine readFile

! ARRAYS

  Subroutine PrintMatrix(xMatrix)
! Prints out a 2D matrix
    Implicit None  ! Force declaration of all variables
! Declare private variables
    Real(kind=DoubleReal), Dimension(:,:) :: xMatrix
    Integer(kind=StandardInteger) :: i, j, k, n
    Character(len=16) :: tempStr
    Character(len=4096) :: printRow
! Print
    Do i=1,size(xMatrix,1)
      Do j=1,4096
        printRow(j:j) = " "
      End Do
      Do j=1,size(xMatrix,2)
        Do k=1,16
          tempStr(k:k) = " "
        End Do
        write(tempStr,"(E10.4)") xMatrix(i,j)
        Do k=1,16
          n = (j-1)*16+k
          printRow(n:n) = tempStr(k:k)
        End Do
      End Do
      print *,trim(printRow)
    End Do
  End Subroutine PrintMatrix

 Subroutine extractArrayColumnDP(inputArray,outputArray,column)
! Extract one column of a 2D dp array array(row,col)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, column
    Real(kind=DoubleReal), Dimension( : , : ), Allocatable :: inputArray
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
! Copy column
    Do i=1,size(inputArray,1)
      outputArray(i) = inputArray(i,column)
    End Do
  End Subroutine extractArrayColumnDP

  Subroutine extractArrayColumnInt(inputArray,outputArray)
! Subroutine extractArrayColumnInt(inputArray,outputArray,column)
! Extract one column of a 2D int array array(row,col)
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger), Dimension( : , : ), Allocatable :: inputArray
    Integer(kind=StandardInteger), Dimension( : ), Allocatable :: outputArray
! Allocate output array
    Allocate(outputArray(1:size(inputArray,1)))
! Copy column
! Do i=1,size(inputArray,1)
!  outputArray(i) = inputArray(i,column)
! End Do
  End Subroutine extractArrayColumnInt

  Subroutine swapArrayRows1D(matrix,rowA,rowB)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
! Set variables
    matH = size(matrix,1)
    matW = 1
! Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
! Allocate arrays
      Allocate(rowAArr(1:matW))
      Allocate(rowBArr(1:matW))
! Swap rows
      Do i=1,matW
        rowAArr(i) = matrix(rowA)
        rowBArr(i) = matrix(rowB)
      End Do
      Do i=1,matW
        matrix(rowA) = rowBArr(i)
        matrix(rowB) = rowAArr(i)
      End Do
    End If
  End Subroutine swapArrayRows1D

  Subroutine swapArrayRows2D(matrix,rowA,rowB)
! Swap rows of square dp matrix
! force declaration of all variables
    Implicit None
! declare private variables
    Integer(kind=StandardInteger) :: i, rowA, rowB, matH, matW
    Real(kind=DoubleReal), Dimension( : , :), Allocatable :: matrix
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowAArr
    Real(kind=DoubleReal), Dimension( : ), Allocatable :: rowBArr
! Set variables
    matH = size(matrix,1)
    matW = size(matrix,2)
! Only do if rows are in the matrix
    If(rowA.ge.1.and.rowA.le.matH.and.rowB.ge.1.and.rowB.le.matH)Then
! Allocate arrays
      Allocate(rowAArr(1:matW))
      Allocate(rowBArr(1:matW))
! Swap rows
      Do i=1,matW
        rowAArr(i) = matrix(rowA,i)
        rowBArr(i) = matrix(rowB,i)
      End Do
      Do i=1,matW
        matrix(rowA,i) = rowBArr(i)
        matrix(rowB,i) = rowAArr(i)
      End Do
    End If
  End Subroutine swapArrayRows2D

  ! Integer, 1D:
  Subroutine swapRows_Integer_1D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:)
    Integer(kind=StandardInteger) :: temp,rowA,rowB
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Integer_1D

! Integer, 2D:
  Subroutine swapRows_Integer_2D(matrix,rowA,rowB)
    Integer(kind=StandardInteger) :: matrix(:,:)
    Integer(kind=StandardInteger) :: temp, i, rowA, rowB
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Integer_2D

! Single, 1D:
  Subroutine swapRows_Single_1D(matrix,rowA,rowB)
    Real(kind=SingleReal) :: matrix(:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=SingleReal) :: temp
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Single_1D

! Single, 2D:
  Subroutine swapRows_Single_2D(matrix,rowA,rowB)
    Real(kind=SingleReal) :: matrix(:,:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=SingleReal) :: temp
    Integer(kind=StandardInteger) :: i
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Single_2D

! Double, 1D:
  Subroutine swapRows_Double_1D(matrix,rowA,rowB)
    Real(kind=DoubleReal) :: matrix(:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=DoubleReal) :: temp
! Swap rows
    temp = matrix(rowA)
    matrix(rowA) = matrix(rowB)
    matrix(rowB) = temp
  End Subroutine swapRows_Double_1D

! Double, 2D:
  Subroutine swapRows_Double_2D(matrix,rowA,rowB)
    Real(kind=DoubleReal) :: matrix(:,:)
    Integer(kind=StandardInteger) :: rowA, rowB
    Real(kind=DoubleReal) :: temp
    Integer(kind=StandardInteger) :: i
! Swap rows
    Do i=1,size(matrix,2)  !Loop through columns
      temp = matrix(rowA,i)
      matrix(rowA,i) = matrix(rowB,i)
      matrix(rowB,i) = temp
    End Do
  End Subroutine swapRows_Double_2D

! Integer, 1D:
  Subroutine sort_Integer_1D(list)
    Integer(kind=StandardInteger) :: list(:)
    Integer(kind=StandardInteger) :: i, sortComplete
! Sort list
    sortComplete = 0
    Do While(sortComplete.eq.0)
      sortComplete = 1
      Do i=2,size(list,1)
        If(list(i-1).gt.list(i))Then
          Call swapRows_Integer_1D(list,i-1,i)
          sortComplete = 0
        End If
      End Do
    End Do
  End Subroutine sort_Integer_1D

! Integer, 2D:
  Subroutine sort_Integer_2D(list, sortRow)
    Integer(kind=StandardInteger) :: list(:,:)
    Integer(kind=StandardInteger) :: i, sortRow, sortComplete
! Sort list
    sortComplete = 0
    Do While(sortComplete.eq.0)
      sortComplete = 1
      Do i=2,size(list,1)
        If(list(i-1,sortRow).gt.list(i,sortRow))Then
          Call swapRows_Integer_2D(list,i-1,i)
          sortComplete = 0
        End If
      End Do
    End Do
  End Subroutine sort_Integer_2D


! TYPE CONVERSION

  Subroutine strToIntArr(stringIn,intArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Integer(kind=StandardInteger), Dimension(:) :: intArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(32) :: intString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    intString = BlankString(intString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        Read(intString,*) intArr(k)
        intString = BlankString(intString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        intString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToIntArr

  Subroutine strToDPArr(stringIn,dpArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Real(kind=DoubleReal), Dimension(:) :: dpArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(32) :: intString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    intString = BlankString(intString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        Read(intString,*) dpArr(k)
        intString = BlankString(intString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        intString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToDPArr

  Subroutine strToStrArr(stringIn,strArr)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Integer(kind=StandardInteger) :: i, j, k
    Character(*), Dimension(:) :: strArr
    Character(*) :: stringIn
    Character(len(stringIn)) :: stringPrep
    Character(Len(strArr)) :: tempString
! Prepare input string
    stringIn = Trim(Adjustl(stringIn))
    stringPrep = BlankString(stringPrep)
! One space only
    j = 0
    Do i=1,len(stringIn)
      If(i.eq.1)Then
        j = j + 1
        stringPrep(j:j) = stringIn(i:i)
      Else
        If(ichar(stringIn(i:i)).eq.32.and.ichar(stringIn(i-1:i-1)).eq.32)Then
! Do not add
        Else
          j = j + 1
          stringPrep(j:j) = stringIn(i:i)
        End If
      End If
    End Do
    tempString = BlankString(tempString)
    j = 0
    k = 0
    Do i=1,len(stringPrep)
      If(ichar(stringPrep(i:i)).eq.32)Then  !Space
        j = 0
        k = k + 1
        strArr(k) = tempString
        tempString = BlankString(tempString)
        If(ichar(stringPrep(i+1:i+1)).eq.32)Then
          Exit
        End If
      Else
        j = j + 1
        tempString(j:j) = stringPrep(i:i)
      End If
    End Do
  End Subroutine strToStrArr

! Time Accumulator subroutine
  Subroutine timeAcc(time,timeStart,timeEnd)
! Take space separated integers and convert to array
    Implicit None   ! Force declaration of all variables
! Declare private variables
    Real(kind=DoubleReal) :: time,timeStart,timeEnd
    time = time + timeEnd - timeStart
  End Subroutine timeAcc


! --------------------------------------------------------------------------------------------------
!    Functions
! --------------------------------------------------------------------------------------------------


  Function FileExists(filePath) Result (boolOut)
    Implicit None   ! Force declaration of all variables
! Private variables
    Character(*) :: filePath
    Logical :: boolOut
! Inquire
    INQUIRE(FILE=trim(filePath), EXIST=boolOut)
  End Function FileExists

  Function CountRowFields(fileRow) Result (fieldCount)
    Implicit None   ! Force declaration of all variables
! In
    Character(*) :: fileRow
! Out
    Integer(kind=StandardInteger) :: fieldCount
! Private variables
    Integer(kind=StandardInteger) :: i
    Logical :: inQuotes
!    Character(Len(fileRow)) :: fileRowTemp

! Blank string
!   fileRowTemp = BlankString(fileRowTemp)
    fieldCount = 1
    inQuotes = .false.
    Do i=1,Len(trim(fileRow))
      If(fileRow(i:i).eq.'"')Then
        If(inQuotes)Then
          inQuotes = .false.
        Else
          inQuotes = .true.
        End If
      End If
      If(.Not.inQuotes)Then
        If(i.gt.1)Then
          If(fileRow(i:i).eq." ".and.fileRow(i-1:i-1).ne." ")Then
            fieldCount = fieldCount + 1
          End If
        End If
      End If
    End Do




  End Function CountRowFields






!---------------------------------------------------------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------------------------------------------------------
End Module general
